// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Babak Farhadi Jahromi, CMC group,
   Ruhr-Universitaet Bochum 

   Based on fix qeq/reaxff by Hasan Metin Aktulga and fix qeq by Ray Shan
------------------------------------------------------------------------- */

#include "fix_qeq_gauss.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "math_const.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixQEqGauss::FixQEqGauss(LAMMPS *lmp, int narg, char **arg) :
  FixQEq(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix qeq/gauss requires atom attribute q");

  ngroup = group->count(igroup);
  if (ngroup == 0) error->all(FLERR,"Fix qeq/gauss group has no atoms");

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  int ntypes = atom->ntypes;
  for (int i = 1; i <= ntypes; i++) {
    if (zeta[i] == 0.0)
      error->all(FLERR,"Invalid param file for fix qeq/gauss");
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  if (atom->nmax > nmax) reallocate_storage();

  if (nlocal > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

  init_matvec();
  matvecs = CG(b_s, s);         // CG on s - parallel
  matvecs += CG(b_t, t);        // CG on t - parallel
  matvecs /= 2;
  calculate_Q();

  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::init_matvec()
{
  compute_H();

  int inum, ii, i;
  int *ilist;

  inum = list->inum;
  ilist = list->ilist;

  for(ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      Hdia_inv[i] = 1. / eta[atom->type[i]];
      b_s[i]      = -(chi[atom->type[i]] + chizj[i]);
      b_t[i]      = -1.0;
      t[i] = t_hist[i][2] + 3 * (t_hist[i][0] - t_hist[i][1]);
      s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
    }
  }

  pack_flag = 2;
  comm->forward_comm_fix(this); //Dist_vector(s);
  pack_flag = 3;
  comm->forward_comm_fix(this); //Dist_vector(t);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::compute_H()
{
  int inum, jnum, *ilist, *jlist, *numneigh, **firstneigh;
  int i, j, ii, jj, itype, jtype;

  double dx, dy, dz, r_sqr, r;
  double zei, zej;

  int *type = atom->type;
  double **x = atom->x;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // fill in the H matrix
  m_fill = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    zei = zeta[itype];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    H.firstnbr[i] = m_fill;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];
      zej = zeta[jtype];
      dx = x[j][0] - x[i][0];
      dy = x[j][1] - x[i][1];
      dz = x[j][2] - x[i][2];
      r_sqr = dx*dx + dy*dy + dz*dz;

      if (r_sqr > cutoff_sq) continue;

      r = sqrt(r_sqr);
      H.jlist[m_fill] = j;
      H.val[m_fill] = calculate_H_dsf(zei, zej, r);
      m_fill++;
    }
    H.numnbrs[i] = m_fill - H.firstnbr[i];
  }

  // DEBUG BABAK
  //for (i=0; i <= 30; ++i)
  //  printf("DEBUG H_mat: %1d %10.5f\n", H.jlist[i], H.val[i]);
  // DEBUG BABAK
  printf("DEBUG: inum: %1d gnum: %1d nlocal: %1d nghost: %1d\n", list->inum, list->gnum, atom->nlocal, atom->nghost);
  if (m_fill >= H.m) {
    char str[128];
    sprintf(str,"H matrix size has been exceeded: m_fill=%d H.m=%d\n",
             m_fill, H.m );
    error->warning(FLERR,str);
    error->all(FLERR,"Fix qeq/gauss has insufficient QEq matrix size");
  }
}

/* ---------------------------------------------------------------------- */

double FixQEqGauss::calculate_H(double zei, double zej, double r)
{
  //double zei2 = zei*zei;
  //double zei2inv = 1.0/zei2;
  //double zej2 = zej*zej;
  //double zej2inv = 1.0/zej2;
  //double zeij = 1.0/sqrt(zei2inv+zej2inv);
  //double erfrze = erf(0.5*r*zeij);
  double sigi = 1.0/(2.0*zei);
  double sigj = 1.0/(2.0*zej);
  double sigij = sqrt(sigi*sigi+sigj*sigj);
  double erfrsiginv = erf(r/sigij);
  double qqrd2e = force->qqrd2e;
  double etmp;
  
  etmp = erfrsiginv/r;
  // DEBUG BABAK
  //printf("DEBUG sigi: %10.5f, sigj: %10.5f, r: %10.5f, etmp: %10.5f, Jij: %10.5f\n", sigi, sigj, r, etmp, 0.5*qqrd2e*etmp);
  
  return 0.5*qqrd2e*etmp;
}

/* ---------------------------------------------------------------------- */

double FixQEqGauss::calculate_H_wolf(double zei, double zej, double r)
{
  double sigi = 1.0/(2.0*zei);
  double sigj = 1.0/(2.0*zej);
  double sigij = sqrt(sigi*sigi+sigj*sigj);
  double siginv = 1.0/sigij;
  double rinv = 1.0/r;
  double rcut = sqrt(cutoff_sq);
  double rcutinv = 1.0/rcut;
  double erfrsiginv = erf(r*siginv);
  double erfrcutsiginv = erf(rcut*siginv);
  double qqrd2e = force->qqrd2e;
  double etmp;
  
  etmp = erfrsiginv*rinv-erfrcutsiginv*rcutinv;
  
  return 0.5*qqrd2e*etmp;
}
/* ---------------------------------------------------------------------- */

double FixQEqGauss::calculate_H_dsf(double zei, double zej, double r)
{
  double sigi = 1.0/(2.0*zei);
  double sigj = 1.0/(2.0*zej);
  double sigij = sqrt(sigi*sigi+sigj*sigj);
  double siginv = 1.0/sigij;
  double rinv = 1.0/r;
  double rcut = sqrt(cutoff_sq);
  double rcutinv = 1.0/rcut;
  double erfrsiginv = erf(r*siginv);
  double erfrcutsiginv = erf(rcut*siginv);
  double preexp = 2.0*siginv/sqrt(MY_PI);
  double expterm = exp(-siginv*siginv*cutoff_sq);
  double qqrd2e = force->qqrd2e;
  double etmp1, etmp2, etmp3;
  
  etmp1 = erfrsiginv*rinv-erfrcutsiginv*rcutinv;
  etmp2 = erfrcutsiginv*rcutinv*rcutinv+preexp*expterm*rcutinv;
  etmp3 = etmp1+etmp2*(r-rcut);
  
  return 0.5*qqrd2e*etmp3;
}

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
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_efield.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "region.h"
#include "respa.h"
#include "text_file_reader.h"
#include "update.h"

#include <cmath>
#include <cstring>
#include <exception>
#include <string>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

static constexpr double EV_TO_KCAL_PER_MOL = 14.4;
static constexpr double SMALL = 1.0e-14;
static constexpr double QSUMSMALL = 0.00001;

/* ---------------------------------------------------------------------- */

FixQEqGauss::FixQEqGauss(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), matvecs(0), pertype_option(nullptr)
{
  scalar_flag = 1;
  extscalar = 0;
  imax = 200;
  maxwarn = 1;

  //if ((narg < 8) || (narg > 12)) error->all(FLERR,"Illegal fix qeq/reaxff command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix qeq/gauss command");

  cutoff = utils::numeric(FLERR,arg[4],false,lmp);
  tolerance = utils::numeric(FLERR,arg[5],false,lmp);
  pertype_option = utils::strdup(arg[6]);

  // dual CG support only available for OPENMP variant
  // check for compatibility is in Fix::post_constructor()

  dual_enabled = 0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dual") == 0) dual_enabled = 1;
    else if (strcmp(arg[iarg],"nowarn") == 0) maxwarn = 0;
    else if (strcmp(arg[iarg],"maxiter") == 0) {
      if (iarg+1 > narg-1)
        error->all(FLERR,"Illegal fix {} command", style);
      imax = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg++;
    } else error->all(FLERR,"Illegal fix {} command", style);
    iarg++;
  }

  nn = n_cap = 0;
  nmax = 0;
  m_fill = m_cap = 0;
  pack_flag = 0;
  s = nullptr;
  t = nullptr;
  nprev = 4;

  Hdia_inv = nullptr;
  b_s = nullptr;
  chi_field = nullptr;
  b_t = nullptr;
  b_prc = nullptr;
  b_prm = nullptr;

  // CG

  p = nullptr;
  q = nullptr;
  r = nullptr;
  d = nullptr;

  // H matrix

  H.firstnbr = nullptr;
  H.numnbrs = nullptr;
  H.jlist = nullptr;
  H.val = nullptr;

  // others
  cutoff_sq = cutoff*cutoff;
  chizj = nullptr;

  // dual CG support
  // Update comm sizes for this fix

  if (dual_enabled) comm_forward = comm_reverse = 2;
  else comm_forward = comm_reverse = 1;

  // perform initial allocation of atom-based arrays
  // register with Atom class

  s_hist = t_hist = nullptr;
  atom->add_callback(Atom::GROW);
}

/* ---------------------------------------------------------------------- */

FixQEqGauss::~FixQEqGauss()
{
  if (copymode) return;

  delete[] pertype_option;

  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,Atom::GROW);

  memory->destroy(s_hist);
  memory->destroy(t_hist);

  FixQEqGauss::deallocate_storage();
  FixQEqGauss::deallocate_matrix();

  memory->destroy(chi);
  memory->destroy(eta);
  memory->destroy(zeta);
  memory->destroy(zcore);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::post_constructor()
{
  grow_arrays(atom->nmax);
  for (int i = 0; i < atom->nmax; i++)
    for (int j = 0; j < nprev; ++j)
      s_hist[i][j] = t_hist[i][j] = 0;

  pertype_parameters(pertype_option);
  if (dual_enabled)
    error->all(FLERR,"Dual keyword only supported with fix qeq/gauss/omp");
}

/* ---------------------------------------------------------------------- */

int FixQEqGauss::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::pertype_parameters(char *arg)
{
  const int ntypes = atom->ntypes;

  memory->create(chi,ntypes+1,"qeq/gauss:chi");
  memory->create(eta,ntypes+1,"qeq/gauss:eta");
  memory->create(zeta,ntypes+1,"qeq/gauss:zeta");
  memory->create(zcore,ntypes+1,"qeq/gauss:zcore");

  if (comm->me == 0) {
    chi[0] = eta[0] = zeta[0] = zcore[0] = 0.0;
    try {
      TextFileReader reader(arg,"qeq/gauss parameter");
      reader.ignore_comments = false;
      for (int i = 1; i <= ntypes; i++) {
        const char *line = reader.next_line();
        if (!line)
          throw TokenizerException("Fix qeq/gauss: Invalid param file format","");
        ValueTokenizer values(line);

        if (values.count() != 5)
          throw TokenizerException("Fix qeq/gauss: Incorrect format of param file","");

        int itype = values.next_int();
        if ((itype < 1) || (itype > ntypes))
          throw TokenizerException("Fix qeq/gauss: invalid atom type in param file",
                                   std::to_string(itype));

        chi[itype] = values.next_double();
        eta[itype] = values.next_double();
        zeta[itype] = values.next_double();
        zcore[itype] = values.next_double();
      }
    } catch (std::exception &e) {
      error->one(FLERR,e.what());
    }
  }

  MPI_Bcast(chi,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(eta,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(zeta,ntypes+1,MPI_DOUBLE,0,world);
  MPI_Bcast(zcore,ntypes+1,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::allocate_storage()
{
  nmax = atom->nmax;

  memory->create(s,nmax,"qeq:s");
  memory->create(t,nmax,"qeq:t");

  memory->create(Hdia_inv,nmax,"qeq:Hdia_inv");
  memory->create(b_s,nmax,"qeq:b_s");
  memory->create(chi_field,nmax,"qeq:chi_field");
  memory->create(b_t,nmax,"qeq:b_t");
  memory->create(b_prc,nmax,"qeq:b_prc");
  memory->create(b_prm,nmax,"qeq:b_prm");

  // dual CG support
  int size = nmax;
  if (dual_enabled) size*= 2;

  memory->create(p,size,"qeq:p");
  memory->create(q,size,"qeq:q");
  memory->create(r,size,"qeq:r");
  memory->create(d,size,"qeq:d");
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::deallocate_storage()
{
  memory->destroy(s);
  memory->destroy(t);

  memory->destroy(Hdia_inv);
  memory->destroy(b_s);
  memory->destroy(b_t);
  memory->destroy(b_prc);
  memory->destroy(b_prm);
  memory->destroy(chi_field);

  memory->destroy(p);
  memory->destroy(q);
  memory->destroy(r);
  memory->destroy(d);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::reallocate_storage()
{
  deallocate_storage();
  allocate_storage();
  init_storage();
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::allocate_matrix()
{
  int i,ii,m;

  int mincap;
  double safezone;

  mincap = MIN_CAP;
  safezone = SAFE_ZONE;

  n_cap = MAX((int)(atom->nlocal * safezone), mincap);

  // determine the total space for the H matrix

  m = 0;
  for (ii = 0; ii < nn; ii++) {
    i = ilist[ii];
    m += numneigh[i];
  }
  m_cap = MAX((int)(m * safezone), mincap * MIN_NBRS);

  H.n = n_cap;
  H.m = m_cap;
  memory->create(H.firstnbr,n_cap,"qeq:H.firstnbr");
  memory->create(H.numnbrs,n_cap,"qeq:H.numnbrs");
  memory->create(H.jlist,m_cap,"qeq:H.jlist");
  memory->create(H.val,m_cap,"qeq:H.val");
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::deallocate_matrix()
{
  memory->destroy(H.firstnbr);
  memory->destroy(H.numnbrs);
  memory->destroy(H.jlist);
  memory->destroy(H.val);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::reallocate_matrix()
{
  deallocate_matrix();
  allocate_matrix();
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::init()
{
  if (!atom->q_flag)
    error->all(FLERR,"Fix {} requires atom attribute q", style);

  if (group->count(igroup) == 0)
    error->all(FLERR,"Fix {} group has no atoms", style);

  // compute net charge and print warning if too large

  double qsum_local = 0.0, qsum = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit)
      qsum_local += atom->q[i];
  }
  MPI_Allreduce(&qsum_local,&qsum,1,MPI_DOUBLE,MPI_SUM,world);

  if ((comm->me == 0) && (fabs(qsum) > QSUMSMALL))
    error->warning(FLERR,"Fix {} group is not charge neutral, net charge = {:.8}", style, qsum);

  // get pointer to fix efield if present. there may be at most one instance of fix efield in use.

  efield = nullptr;
  auto fixes = modify->get_fix_by_style("^efield");
  if (fixes.size() == 1) efield = dynamic_cast<FixEfield *>( fixes.front());
  else if (fixes.size() > 1)
    error->all(FLERR, "There may be only one fix efield instance used with fix {}", style);

  // ensure that fix efield is properly initialized before accessing its data and check some settings
  if (efield) {
    efield->init();
    if (strcmp(update->unit_style,"real") != 0)
      error->all(FLERR,"Must use unit_style real with fix {} and external fields", style);
    if (efield->varflag != FixEfield::CONSTANT)
      error->all(FLERR,"Cannot (yet) use fix {} with variable efield", style);

    if (((fabs(efield->ex) > SMALL) && domain->xperiodic) ||
         ((fabs(efield->ey) > SMALL) && domain->yperiodic) ||
         ((fabs(efield->ez) > SMALL) && domain->zperiodic))
      error->all(FLERR,"Must not have electric field component in direction of periodic "
                       "boundary when using charge equilibration with ReaxFF.");
  }

  // we need a half neighbor list w/ Newton off
  // built whenever re-neighboring occurs

  neighbor->add_request(this, NeighConst::REQ_NEWTON_OFF);

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = (dynamic_cast<Respa *>( update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

double FixQEqGauss::compute_scalar()
{
  return matvecs/2.0;
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::setup_pre_force(int vflag)
{
  nn = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  deallocate_storage();
  allocate_storage();

  init_storage();

  deallocate_matrix();
  allocate_matrix();

  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel < nlevels_respa-1) return;
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::min_setup_pre_force(int vflag)
{
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::init_storage()
{
  if (efield) get_chi_field();

  for (int ii = 0; ii < nn; ii++) {
    int i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      Hdia_inv[i] = 1. / eta[atom->type[i]];
      b_s[i] = -chi[atom->type[i]];
      if (efield) b_s[i] -= chi_field[i];
      b_t[i] = -1.0;
      b_prc[i] = 0;
      b_prm[i] = 0;
      s[i] = t[i] = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  int n = atom->nlocal;

  nn = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // grow arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) reallocate_storage();
  if (n > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

  if (efield) get_chi_field();

  init_matvec();

  matvecs_s = CG(b_s, s);       // CG on s - parallel
  matvecs_t = CG(b_t, t);       // CG on t - parallel
  matvecs = matvecs_s + matvecs_t;

  calculate_Q();
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::pre_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::init_matvec()
{
  /* fill-in H matrix */
  compute_H();

  int ii, i;

  for (ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {

      /* init pre-conditioner for H and init solution vectors */
      Hdia_inv[i] = 1. / eta[atom->type[i]];
      b_s[i]      = -chi[atom->type[i]];
      if (efield) b_s[i] -= chi_field[i];
      b_t[i]      = -1.0;

      /* quadratic extrapolation for s & t from previous solutions */
      t[i] = t_hist[i][2] + 3 * (t_hist[i][0] - t_hist[i][1]);

      /* cubic extrapolation for s & t from previous solutions */
      s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
    }
  }

  pack_flag = 2;
  comm->forward_comm(this); //Dist_vector(s);
  pack_flag = 3;
  comm->forward_comm(this); //Dist_vector(t);
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::compute_H()
{
  int jnum;
  int i, j, ii, jj, itype, jtype, flag;
  double dx, dy, dz, r_sqr, zei, zej;
  constexpr double EPSILON = 0.0001;

  int *type = atom->type;
  tagint *tag = atom->tag;
  double **x = atom->x;
  int *mask = atom->mask;

  // fill in the H matrix
  m_fill = 0;
  r_sqr = 0;
  for (ii = 0; ii < nn; ii++) {
    i = ilist[ii];
    itype = type[i];
    zei = zeta[itype];
    if (mask[i] & groupbit) {
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
        r_sqr = SQR(dx) + SQR(dy) + SQR(dz);

        flag = 0;
        // DEBUG BABAK
        //printf("DEBUG i: %1d, j: %1d, r_sqr: %10.5f\n",i, j, r_sqr);
        if (r_sqr <= cutoff_sq) {
          if (j < atom->nlocal) flag = 1;
          else if (tag[i] < tag[j]) flag = 1;
          else if (tag[i] == tag[j]) {
            if (dz > EPSILON) flag = 1;
            else if (fabs(dz) < EPSILON) {
              if (dy > EPSILON) flag = 1;
              else if (fabs(dy) < EPSILON && dx > EPSILON)
                flag = 1;
            }
          }
        }

        if (flag) {
          H.jlist[m_fill] = j;
          H.val[m_fill] =  calculate_H_dsf(zei, zej, sqrt(r_sqr));
          m_fill++;
        }
      }
      H.numnbrs[i] = m_fill - H.firstnbr[i];
    }
  }

  if (m_fill >= H.m)
    error->all(FLERR,fmt::format("Fix qeq/reaxff H matrix size has been "
                                 "exceeded: m_fill={} H.m={}\n", m_fill, H.m));
}

/* ---------------------------------------------------------------------- */

double FixQEqGauss::calculate_H(double zei, double zej, double r)
{
  double sigi = 1.0/(2.0*zei);
  double sigj = 1.0/(2.0*zej);
  double sigij = sqrt(sigi*sigi+sigj*sigj);
  double erfrsiginv = erf(r/sigij);
  double qqrd2e = force->qqrd2e;
  double etmp;
  
  etmp = erfrsiginv/r;
  
  return qqrd2e*etmp;
}

/* ---------------------------------------------------------------------- */

double FixQEqGauss::calculate_H_wolf(double zei, double zej, double r)
{
  double sigi = 1.0/(2.0*zei);
  double sigj = 1.0/(2.0*zej);
  double sigij = sqrt(sigi*sigi+sigj*sigj);
  double siginv = 1.0/sigij;
  double rinv = 1.0/r;
  double rcut = cutoff;
  double rcutinv = 1.0/rcut;
  double erfrsiginv = erf(r*siginv);
  double erfrcutsiginv = erf(rcut*siginv);
  double qqrd2e = force->qqrd2e;
  double etmp;
  
  etmp = erfrsiginv*rinv-erfrcutsiginv*rcutinv;
  
  return qqrd2e*etmp;
}
/* ---------------------------------------------------------------------- */

double FixQEqGauss::calculate_H_dsf(double zei, double zej, double r)
{
  double sigi = 1.0/(2.0*zei);
  double sigj = 1.0/(2.0*zej);
  double sigij = sqrt(sigi*sigi+sigj*sigj);
  double siginv = 1.0/sigij;
  double rinv = 1.0/r;
  double rcut = cutoff;
  double rcutinv = 1.0/rcut;
  double erfrsiginv = erf(r*siginv);
  double erfrcutsiginv = erf(rcut*siginv);
  double preexp = 2.0*siginv/sqrt(MY_PI);
  double expterm = exp(-siginv*siginv*cutoff_sq);
  double qqrd2e = force->qqrd2e;
  double etmp1, etmp2, etmp3;
  
  etmp1 = erfrsiginv*rinv-erfrcutsiginv*rcutinv;
  etmp2 = erfrcutsiginv*rcutinv*rcutinv-preexp*expterm*rcutinv;
  etmp3 = etmp1+etmp2*(r-rcut);
  // DEBUG BABAK
  //printf("DEBUG sigi: %10.5f, sigj: %10.5f, r: %10.5f, etmp: %10.5f, Jij: %10.5f\n", sigi, sigj, r, etmp3, qqrd2e*etmp3);
  
  return qqrd2e*etmp3;
}

/* ---------------------------------------------------------------------- */

int FixQEqGauss::CG(double *b, double *x)
{
  int  i, j;
  double tmp, alpha, beta, b_norm;
  double sig_old, sig_new;

  int jj;

  pack_flag = 1;
  sparse_matvec(&H, x, q);
  comm->reverse_comm(this); //Coll_Vector(q);

  vector_sum(r , 1.,  b, -1., q, nn);

  for (jj = 0; jj < nn; ++jj) {
    j = ilist[jj];
    if (atom->mask[j] & groupbit)
      d[j] = r[j] * Hdia_inv[j]; //pre-condition
  }

  b_norm = parallel_norm(b, nn);
  sig_new = parallel_dot(r, d, nn);

  for (i = 1; i < imax && sqrt(sig_new) / b_norm > tolerance; ++i) {
    comm->forward_comm(this); //Dist_vector(d);
    sparse_matvec(&H, d, q);
    comm->reverse_comm(this); //Coll_vector(q);

    tmp = parallel_dot(d, q, nn);
    alpha = sig_new / tmp;

    vector_add(x, alpha, d, nn);
    vector_add(r, -alpha, q, nn);

    // pre-conditioning
    for (jj = 0; jj < nn; ++jj) {
      j = ilist[jj];
      if (atom->mask[j] & groupbit)
        p[j] = r[j] * Hdia_inv[j];
    }

    sig_old = sig_new;
    sig_new = parallel_dot(r, p, nn);

    beta = sig_new / sig_old;
    vector_sum(d, 1., p, beta, d, nn);
  }

  if ((i >= imax) && maxwarn && (comm->me == 0))
    error->warning(FLERR,fmt::format("Fix qeq/gauss CG convergence failed "
                                     "after {} iterations at step {}",
                                     i,update->ntimestep));
  return i;
}


/* ---------------------------------------------------------------------- */

void FixQEqGauss::sparse_matvec(sparse_matrix *A, double *x, double *b)
{
  int i, j, itr_j;
  int ii;

  for (ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      b[i] = eta[atom->type[i]] * x[i];
  }

  int nall = atom->nlocal + atom->nghost;
  for (i = atom->nlocal; i < nall; ++i)
      b[i] = 0;

  for (ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      for (itr_j=A->firstnbr[i]; itr_j<A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
        j = A->jlist[itr_j];
        b[i] += A->val[itr_j] * x[j];
        b[j] += A->val[itr_j] * x[i];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::calculate_Q()
{
  int i, k;
  double u, s_sum, t_sum;
  double *q = atom->q;

  int ii;

  s_sum = parallel_vector_acc(s, nn);
  t_sum = parallel_vector_acc(t, nn);
  u = s_sum / t_sum;

  for (ii = 0; ii < nn; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      q[i] = s[i] - u * t[i];

      /* backup s & t */
      for (k = nprev-1; k > 0; --k) {
        s_hist[i][k] = s_hist[i][k-1];
        t_hist[i][k] = t_hist[i][k-1];
      }
      s_hist[i][0] = s[i];
      t_hist[i][0] = t[i];
    }
  }

  pack_flag = 4;
  comm->forward_comm(this); //Dist_vector(atom->q);
}

/* ---------------------------------------------------------------------- */

int FixQEqGauss::pack_forward_comm(int n, int *list, double *buf,
                                  int /*pbc_flag*/, int * /*pbc*/)
{
  int m;

  if (pack_flag == 1)
    for (m = 0; m < n; m++) buf[m] = d[list[m]];
  else if (pack_flag == 2)
    for (m = 0; m < n; m++) buf[m] = s[list[m]];
  else if (pack_flag == 3)
    for (m = 0; m < n; m++) buf[m] = t[list[m]];
  else if (pack_flag == 4)
    for (m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  else if (pack_flag == 5) {
    m = 0;
    for (int i = 0; i < n; i++) {
      int j = 2 * list[i];
      buf[m++] = d[j];
      buf[m++] = d[j+1];
    }
    return m;
  }
  return n;
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

  if (pack_flag == 1)
    for (m = 0, i = first; m < n; m++, i++) d[i] = buf[m];
  else if (pack_flag == 2)
    for (m = 0, i = first; m < n; m++, i++) s[i] = buf[m];
  else if (pack_flag == 3)
    for (m = 0, i = first; m < n; m++, i++) t[i] = buf[m];
  else if (pack_flag == 4)
    for (m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
  else if (pack_flag == 5) {
    int last = first + n;
    m = 0;
    for (i = first; i < last; i++) {
      int j = 2 * i;
      d[j] = buf[m++];
      d[j+1] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixQEqGauss::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  if (pack_flag == 5) {
    m = 0;
    int last = first + n;
    for (i = first; i < last; i++) {
      int indxI = 2 * i;
      buf[m++] = q[indxI];
      buf[m++] = q[indxI+1];
    }
    return m;
  } else {
    for (m = 0, i = first; m < n; m++, i++) buf[m] = q[i];
    return n;
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::unpack_reverse_comm(int n, int *list, double *buf)
{
  if (pack_flag == 5) {
    int m = 0;
    for (int i = 0; i < n; i++) {
      int indxI = 2 * list[i];
      q[indxI] += buf[m++];
      q[indxI+1] += buf[m++];
    }
  } else {
    for (int m = 0; m < n; m++) q[list[m]] += buf[m];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixQEqGauss::memory_usage()
{
  double bytes;

  bytes = (double)atom->nmax*nprev*2 * sizeof(double); // s_hist & t_hist
  bytes += (double)atom->nmax*11 * sizeof(double); // storage
  bytes += (double)n_cap*2 * sizeof(int); // matrix...
  bytes += (double)m_cap * sizeof(int);
  bytes += (double)m_cap * sizeof(double);

  if (dual_enabled)
    bytes += (double)atom->nmax*4 * sizeof(double); // double size for q, d, r, and p

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate fictitious charge arrays
------------------------------------------------------------------------- */

void FixQEqGauss::grow_arrays(int nmax)
{
  memory->grow(s_hist,nmax,nprev,"qeq:s_hist");
  memory->grow(t_hist,nmax,nprev,"qeq:t_hist");
}

/* ----------------------------------------------------------------------
   copy values within fictitious charge arrays
------------------------------------------------------------------------- */

void FixQEqGauss::copy_arrays(int i, int j, int /*delflag*/)
{
  for (int m = 0; m < nprev; m++) {
    s_hist[j][m] = s_hist[i][m];
    t_hist[j][m] = t_hist[i][m];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixQEqGauss::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nprev; m++) buf[m] = s_hist[i][m];
  for (int m = 0; m < nprev; m++) buf[nprev+m] = t_hist[i][m];
  return nprev*2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixQEqGauss::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nprev; m++) s_hist[nlocal][m] = buf[m];
  for (int m = 0; m < nprev; m++) t_hist[nlocal][m] = buf[nprev+m];
  return nprev*2;
}

/* ---------------------------------------------------------------------- */

double FixQEqGauss::parallel_norm(double *v, int n)
{
  int  i;
  double my_sum, norm_sqr;

  int ii;

  my_sum = 0.0;
  norm_sqr = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_sum += SQR(v[i]);
  }

  MPI_Allreduce(&my_sum, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world);

  return sqrt(norm_sqr);
}

/* ---------------------------------------------------------------------- */

double FixQEqGauss::parallel_dot(double *v1, double *v2, int n)
{
  int  i;
  double my_dot, res;

  int ii;

  my_dot = 0.0;
  res = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_dot += v1[i] * v2[i];
  }

  MPI_Allreduce(&my_dot, &res, 1, MPI_DOUBLE, MPI_SUM, world);

  return res;
}

/* ---------------------------------------------------------------------- */

double FixQEqGauss::parallel_vector_acc(double *v, int n)
{
  int  i;
  double my_acc, res;

  int ii;

  my_acc = 0.0;
  res = 0.0;
  for (ii = 0; ii < n; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit)
      my_acc += v[i];
  }

  MPI_Allreduce(&my_acc, &res, 1, MPI_DOUBLE, MPI_SUM, world);

  return res;
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::vector_sum(double* dest, double c, double* v,
                                double d, double* y, int k)
{
  int kk;

  for (--k; k>=0; --k) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] = c * v[kk] + d * y[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::vector_add(double* dest, double c, double* v, int k)
{
  int kk;

  for (--k; k>=0; --k) {
    kk = ilist[k];
    if (atom->mask[kk] & groupbit)
      dest[kk] += c * v[kk];
  }
}

/* ---------------------------------------------------------------------- */

void FixQEqGauss::get_chi_field()
{
  memset(&chi_field[0],0,atom->nmax*sizeof(double));
  if (!efield) return;

  const auto x = (const double * const *)atom->x;
  const int *mask = atom->mask;
  const imageint *image = atom->image;
  const int nlocal = atom->nlocal;


  // update electric field region if necessary

  Region *region = efield->region;
  if (region) region->prematch();

  // efield energy is in real units of kcal/mol/angstrom, need to convert to eV

  const double factor = -1.0/force->qe2f;

  // currently we only support constant efield
  // atom selection is for the group of fix efield

  if (efield->varflag == FixEfield::CONSTANT) {
    double unwrap[3];
    const double fx = efield->ex;
    const double fy = efield->ey;
    const double fz = efield->ez;
    const int efgroupbit = efield->groupbit;

    // charge interactions
    // force = qE, potential energy = F dot x in unwrapped coords

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & efgroupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);
        chi_field[i] = factor*(fx*unwrap[0] + fy*unwrap[1] + fz*unwrap[2]);
      }
    }
  }
}
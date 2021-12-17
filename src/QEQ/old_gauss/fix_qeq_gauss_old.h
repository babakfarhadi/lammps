/* -*- c++ -*- ----------------------------------------------------------
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

   Therefore, please cite the related publication:
   H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
   "Parallel Reactive Molecular Dynamics: Numerical Methods and
   Algorithmic Techniques", Parallel Computing, in press.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(qeq/gauss,FixQEqGauss);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_GAUSS_H
#define LMP_FIX_QEQ_GAUSS_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqGauss : public FixQEq {
 public:
  FixQEqGauss(class LAMMPS *, int, char **);
  ~FixQEqGauss() {}
  void init();
  void pre_force(int);

 private:
  void init_matvec();
  void compute_H();
  double calculate_H(double, double, double);
  double calculate_H_wolf(double, double, double);
  double calculate_H_dsf(double, double, double);

};
}
#endif
#endif

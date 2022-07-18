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

#ifdef FIX_CLASS
// clang-format off
FixStyle(acks2/gauss,FixACKS2Gauss);
// clang-format on
#else

#ifndef LMP_FIX_ACKS2_GAUSS_H
#define LMP_FIX_ACKS2_GAUSS_H

#include "fix_qeq_gauss.h"

namespace LAMMPS_NS {

class FixACKS2Gauss : public FixQEqGauss {
 public:
  FixACKS2Gauss(class LAMMPS *, int, char **);
  ~FixACKS2Gauss() override;
  void post_constructor() override;
  void init() override;
  void init_storage() override;
  void pre_force(int) override;

  double get_cutoff() { return cutoff; }
  double *get_s() { return s; }
  double *get_u() { return u; }
  double *get_chi() { return chi; }
  double *get_eta() { return eta; }
  double *get_Xij() { return Xij; }
  double *get_X_diag() { return X_diag; }

 protected:
  int NN, last_rows_rank, last_rows_flag;

  double **s_hist_X, **s_hist_last;
  double *Xij;    // acks2 parameters

  sparse_matrix X;
  double *Xdia_inv;
  double *X_diag;
  double *u;

  //BiCGStab storage
  double *g, *q_hat, *r_hat, *y, *z;

  void pertype_parameters(char *) override;
  void init_bondcut();
  void allocate_storage() override;
  void deallocate_storage() override;
  void allocate_matrix() override;
  void deallocate_matrix() override;

  void init_matvec() override;
  void compute_X();
  double calculate_X_bonded(double, double, double);
  double calculate_X_nonbonded(double, double, double);
  void calculate_Q() override;

  int BiCGStab(double *, double *);
  void sparse_matvec_acks2(sparse_matrix *, sparse_matrix *, double *, double *);

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void more_forward_comm(double *);
  void more_reverse_comm(double *);
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  double parallel_norm(double *, int) override;
  double parallel_dot(double *, double *, int) override;
  double parallel_vector_acc(double *, int) override;

  void vector_sum(double *, double, double *, double, double *, int) override;
  void vector_add(double *, double, double *, int) override;
  void vector_copy(double *, double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

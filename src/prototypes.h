#ifndef PROTOTYPES_H
#define PROTOTYPES_H 

#include "athena.h"

void add1_2name(char file[9]);

void binary_dump(struct grid_block *agrid);

void data_output(struct grid_block *agrid, int hst_flag, int hdf_flag,
   int bin_flag);

void Eigensystem4_isothermal_hydro_InPrimVars(Real d, Real v1, Real v2, Real v3,
  Real eigenvalues[NVAR],
  Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_adiabatic_hydro_InPrimVars(Real d, Real v1, Real v2, Real v3,
  Real p, Real eigenvalues[NVAR],
  Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_isothermal_mhd_InPrimVars(Real d, Real v1, Real v2, Real v3,
  Real b1, Real b2, Real b3, Real eigenvalues[NVAR],
  Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_adiabatic_mhd_InPrimVars(Real d, Real v1, Real v2, Real v3,
  Real p, Real b1, Real b2, Real b3, Real eigenvalues[NVAR],
  Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_isothermal_hydro_InConsVars(Real v1, Real v2, Real v3,
  Real eigenvalues[NVAR],
  Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_adiabatic_hydro_InConsVars(Real v1, Real v2, Real v3, Real h,
  Real eigenvalues[NVAR],
  Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_isothermal_mhd_InConsVars(Real d, Real v1, Real v2, Real v3,
  Real b1, Real b2, Real b3, Real x, Real y, Real eigenvalues[NVAR],
  Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_adiabatic_mhd_InConsVars(Real d, Real v1, Real v2, Real v3,
  Real h, Real b1, Real b2, Real b3, Real x, Real y,
  Real eigenvalues[NVAR],
  Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR]);

Real init_dt(struct grid_block *agrid);

void init_grid_block(FILE *p_input_file, FILE *p_output_file,
   struct grid_block *agrid);

void integrate_1step(struct grid_block *agrid);

Real lr_states(Real u[NXMAX][NVAR], Real b_parallel[NXMAX+1], Real dtodx,
   int ibegin, int iend, Real wl[NXMAX][NVAR], Real wr[NXMAX][NVAR]);

void printd(struct grid_block *agrid);

Real roe_fluxes(Real wl[NXMAX][NVAR], Real wr[NXMAX][NVAR],
  Real b_parallel[NXMAX+1],int ibegin,int iend, Real f[NXMAX][NVAR]);

void set_bval_arrays(struct grid_block *agrid, struct bval_array *bval);

void set_ghost_zones(struct grid_block *agrid, struct bval_array *bval);

#endif /* PROTOTYPES_H */

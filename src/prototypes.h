#ifndef PROTOTYPES_H
#define PROTOTYPES_H 

void add1_2name(char file[9]);

void binary_dump(struct grid_block *agrid);

void data_output(struct grid_block *agrid, int hst_flag, int hdf_flag,
   int bin_flag);

void Eigensystem4_isothermal_hydro_InPrimVars(REAL d, REAL v1, REAL v2, REAL v3,
  REAL eigenvalues[NVAR],
  REAL right_eigenmatrix[NVAR][NVAR], REAL left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_adiabatic_hydro_InPrimVars(REAL d, REAL v1, REAL v2, REAL v3,
  REAL p, REAL eigenvalues[NVAR],
  REAL right_eigenmatrix[NVAR][NVAR], REAL left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_isothermal_mhd_InPrimVars(REAL d, REAL v1, REAL v2, REAL v3,
  REAL b1, REAL b2, REAL b3, REAL eigenvalues[NVAR],
  REAL right_eigenmatrix[NVAR][NVAR], REAL left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_adiabatic_mhd_InPrimVars(REAL d, REAL v1, REAL v2, REAL v3,
  REAL p, REAL b1, REAL b2, REAL b3, REAL eigenvalues[NVAR],
  REAL right_eigenmatrix[NVAR][NVAR], REAL left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_isothermal_hydro_InConsVars(REAL v1, REAL v2, REAL v3,
  REAL eigenvalues[NVAR],
  REAL right_eigenmatrix[NVAR][NVAR], REAL left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_adiabatic_hydro_InConsVars(REAL v1, REAL v2, REAL v3, REAL h,
  REAL eigenvalues[NVAR],
  REAL right_eigenmatrix[NVAR][NVAR], REAL left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_isothermal_mhd_InConsVars(REAL d, REAL v1, REAL v2, REAL v3,
  REAL b1, REAL b2, REAL b3, REAL x, REAL y, REAL eigenvalues[NVAR],
  REAL right_eigenmatrix[NVAR][NVAR], REAL left_eigenmatrix[NVAR][NVAR]);

void Eigensystem4_adiabatic_mhd_InConsVars(REAL d, REAL v1, REAL v2, REAL v3,
  REAL h, REAL b1, REAL b2, REAL b3, REAL x, REAL y,
  REAL eigenvalues[NVAR],
  REAL right_eigenmatrix[NVAR][NVAR], REAL left_eigenmatrix[NVAR][NVAR]);

REAL init_dt(struct grid_block *agrid);

void init_grid_block(FILE *p_input_file, FILE *p_output_file,
   struct grid_block *agrid);

void integrate_1step(struct grid_block *agrid);

REAL lr_states(REAL u[NXMAX][NVAR], REAL b_parallel[NXMAX+1], REAL dtodx,
   int ibegin, int iend, REAL wl[NXMAX][NVAR], REAL wr[NXMAX][NVAR]);

void printd(struct grid_block *agrid);

REAL roe_fluxes(REAL wl[NXMAX][NVAR], REAL wr[NXMAX][NVAR],
  REAL b_parallel[NXMAX+1],int ibegin,int iend, REAL f[NXMAX][NVAR]);

void set_bval_arrays(struct grid_block *agrid, struct bval_array *bval);

void set_ghost_zones(struct grid_block *agrid, struct bval_array *bval);

#endif /* PROTOTYPES_H */

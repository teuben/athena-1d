#ifndef ATHENA_H
#define ATHENA_H 

#include "athena.def"

typedef double Real;

struct bval_array{
  Real uiib[4][NVAR], uoib[4][NVAR], bxiib[4], bxoib[4];
};

struct grid_block{
  Real u[NX1][NVAR],bx[NX1+1];
  Real dx,dt,time;
  Real t_bin,dt_bin,t_hst,dt_hst,t_hdf,dt_hdf;
  struct bval_array boundary_values;
  int is,ie;
  int niib,noib;
  int ncycles;
  char bin_file[9],hst_file[9],hdf_file[9];
};

#endif /* ATHENA_H */

#ifndef ATHENA_H
#define ATHENA_H 

#include "athena.def"

typedef double Real;

typedef struct bval_array{
#if defined ONE_D
  Real uiib[4][NVAR], uoib[4][NVAR];
#ifdef MHD
  Real bxiib[4], bxoib[4];
#endif /* MHD */

#elif defined TWO_D
#error: 2D bval_array is not yet defined

#else /* THREE_D */
#error: 3D bval_array is not yet defined

#endif
}Bval_Array;



typedef struct grid_block{
#if defined ONE_D
  Real dx;
  Real u[NX1][NVAR];
#ifdef MHD
  Real bx[NX1];
#endif /* MHD */

#elif defined TWO_D
  Real dx,dy;
  Real u[NX1][NX2][NVAR];
#ifdef MHD
  Real bx[NX1][NX2],by[NX1][NX2];
#endif /* MHD */

#else /* THREE_D */
  Real dx,dy,dz;
  Real u[NX1][NX2][NX3][NVAR];
#ifdef MHD
  Real bx[NX1][NX2][NX3],by[NX1][NX2][NX3],bz[NX1][NX2][NX3];
#endif /* MHD */
#endif
  Real dt,time;
  Real t_bin,dt_bin,t_hst,dt_hst,t_hdf,dt_hdf;
  struct bval_array boundary_values;
  int is,ie;
  int niib,noib;
  int ncycles;
  char bin_file[9],hst_file[9],hdf_file[9];
}Grid_Block;

#endif /* ATHENA_H */

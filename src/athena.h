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
  Real uiib1[4][NX2][NVAR], uoib1[4][NX2][NVAR];
  Real uiib2[NX1][4][NVAR], uoib2[NX1][4][NVAR];
#ifdef MHD
  Real bxiib1[4][NX2], bxoib1[4][NX2];
  Real byiib1[4][NX2], byoib1[4][NX2];

  Real bxiib2[NX1][4], bxoib2[NX1][4];
  Real byiib2[NX1][4], byoib2[NX1][4];
#endif /* MHD */

#else /* THREE_D */
  Real uiib1[4][NX2][NX3][NVAR], uoib1[4][NX2][NX3][NVAR];
  Real uiib2[NX1][4][NX3][NVAR], uoib2[NX1][4][NX3][NVAR];
  Real uiib3[NX1][NX2][4][NVAR], uoib3[NX1][NX2][4][NVAR];
#ifdef MHD
  Real bxiib1[4][NX2][NX3], bxoib1[4][NX2][NX3];
  Real byiib1[4][NX2][NX3], byoib1[4][NX2][NX3];
  Real bziib1[4][NX2][NX3], bzoib1[4][NX2][NX3];

  Real bxiib2[NX1][4][NX3], bxoib2[NX1][4][NX3];
  Real byiib2[NX1][4][NX3], byoib2[NX1][4][NX3];
  Real bziib2[NX1][4][NX3], bzoib2[NX1][4][NX3];

  Real bxiib3[NX1][NX2][4], bxoib3[NX1][NX2][4];
  Real byiib3[NX1][NX2][4], byoib3[NX1][NX2][4];
  Real bziib3[NX1][NX2][4], bzoib3[NX1][NX2][4];
#endif /* MHD */

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

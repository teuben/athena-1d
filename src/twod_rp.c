/*============================================================================*/
/*///////////////////////////// Function twod_rp \\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
/*                                                                            */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

/* Initial States Positions:
   y
   |
   ---------
   | 2 | 1 |
   ---------
   | 3 | 4 |
   ----------- x
*/


/* Problem generator for 2-D Riemann problems */
void twod_rp(FILE *p_input_file, struct grid_block *agrid){

  char buf[120];
  int i,j;
  int is1,ie1,is2,ie2;
  int im1,im2;
  Real d1,p1,u1,v1,w1,d2,p2,u2,v2,w2;
  Real d3,p3,u3,v3,w3,d4,p4,u4,v4,w4;
#ifdef MHD
  Real bx1,by1,bz1,bx2,by2,bz2;
  Real bx3,by3,bz3,bx4,by4,bz4;
#endif /* MHD */
  /*==========================================================================*/

  /* Read d1,p1,u1,v1,w1,bx1,by1,bz1 */

  if (fscanf(p_input_file,"%lg %[^\n]\n", &d1, buf) != 2) {
    printf("Error reading d1\n"); exit(EXIT_FAILURE);
  }
#ifdef ADIABATIC
  if (fscanf(p_input_file,"%lg %[^\n]\n", &p1, buf) != 2) {
    printf("Error reading p1\n"); exit(EXIT_FAILURE);
  }
#endif
  if (fscanf(p_input_file,"%lg %[^\n]\n", &u1, buf) != 2) {
    printf("Error reading u1\n"); exit(EXIT_FAILURE);
   }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &v1, buf) != 2) {
    printf("Error reading v1\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &w1, buf) != 2) {
    printf("Error reading w1\n"); exit(EXIT_FAILURE);
  }
#ifdef MHD
  if (fscanf(p_input_file,"%lg %[^\n]\n", &bx1, buf) != 2) {
    printf("Error reading bx1\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &by1, buf) != 2) {
    printf("Error reading by1\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &bz1, buf) != 2) {
    printf("Error reading bz1\n"); exit(EXIT_FAILURE);
  }
#endif

  /* Read d2,p2,u2,v2,w2,bx2,by2,bz2 */

  if (fscanf(p_input_file,"%lg %[^\n]\n", &d2, buf) != 2) {
    printf("Error reading d2\n"); exit(EXIT_FAILURE);
  }
#ifdef ADIABATIC
  if (fscanf(p_input_file,"%lg %[^\n]\n", &p2, buf) != 2) {
    printf("Error reading p2\n"); exit(EXIT_FAILURE);
  }
#endif
  if (fscanf(p_input_file,"%lg %[^\n]\n", &u2, buf) != 2) {
    printf("Error reading u2\n"); exit(EXIT_FAILURE);
   }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &v2, buf) != 2) {
    printf("Error reading v2\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &w2, buf) != 2) {
    printf("Error reading w2\n"); exit(EXIT_FAILURE);
  }
#ifdef MHD
  if (fscanf(p_input_file,"%lg %[^\n]\n", &bx2, buf) != 2) {
    printf("Error reading bx2\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &by2, buf) != 2) {
    printf("Error reading by2\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &bz2, buf) != 2) {
    printf("Error reading bz2\n"); exit(EXIT_FAILURE);
  }
#endif

  /* Read d3,p3,u3,v3,w3,bx3,by3,bz3 */

  if (fscanf(p_input_file,"%lg %[^\n]\n", &d3, buf) != 2) {
    printf("Error reading d3\n"); exit(EXIT_FAILURE);
  }
#ifdef ADIABATIC
  if (fscanf(p_input_file,"%lg %[^\n]\n", &p3, buf) != 2) {
    printf("Error reading p3\n"); exit(EXIT_FAILURE);
  }
#endif
  if (fscanf(p_input_file,"%lg %[^\n]\n", &u3, buf) != 2) {
    printf("Error reading u3\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &v3, buf) != 2) {
    printf("Error reading v3\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &w3, buf) != 2) {
    printf("Error reading w3\n"); exit(EXIT_FAILURE);
  }
#ifdef MHD
  if (fscanf(p_input_file,"%lg %[^\n]\n", &bx3, buf) != 2) {
    printf("Error reading bx3\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &by3, buf) != 2) {
    printf("Error reading by3\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &bz3, buf) != 2) {
    printf("Error reading bz3\n"); exit(EXIT_FAILURE);
  }
#endif

  /* Read d4,p4,u4,v4,w4,bx4,by4,bz4 */

  if (fscanf(p_input_file,"%lg %[^\n]\n", &d4, buf) != 2) {
    printf("Error reading d4\n"); exit(EXIT_FAILURE);
  }
#ifdef ADIABATIC
  if (fscanf(p_input_file,"%lg %[^\n]\n", &p4, buf) != 2) {
    printf("Error reading p4\n"); exit(EXIT_FAILURE);
  }
#endif
  if (fscanf(p_input_file,"%lg %[^\n]\n", &u4, buf) != 2) {
    printf("Error reading u4\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &v4, buf) != 2) {
    printf("Error reading v4\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &w4, buf) != 2) {
    printf("Error reading w4\n"); exit(EXIT_FAILURE);
  }
#ifdef MHD
  if (fscanf(p_input_file,"%lg %[^\n]\n", &bx4, buf) != 2) {
    printf("Error reading bx4\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &by4, buf) != 2) {
    printf("Error reading by4\n"); exit(EXIT_FAILURE);
  }
  if (fscanf(p_input_file,"%lg %[^\n]\n", &bz4, buf) != 2) {
    printf("Error reading bz4\n"); exit(EXIT_FAILURE);
  }
#endif


  /* setup dependent variables for shocktube in X-direction */

  is1 = agrid->is1;  ie1 = agrid->ie1;
  im1 = is1+((ie1-is1+1)/2);
  is2 = agrid->is2;  ie2 = agrid->ie2;
  im2 = is2+((ie2-is2+1)/2);

  for(i=im1; i<=ie1+2; i++){
    for(j=im2; j<=ie2+2; j++){
      agrid->u[i][j][0] = d1;
      agrid->u[i][j][1] = u1*d1;
      agrid->u[i][j][2] = v1*d1;
      agrid->u[i][j][3] = w1*d1;
#ifdef MHD
      agrid->bx[i][j] = bx1;
      agrid->by[i][j] = by1;
      agrid->u[i][j][NVAR-3] = bx1;
      agrid->u[i][j][NVAR-2] = by1;
      agrid->u[i][j][NVAR-1] = bz1;
#endif
#ifdef ADIABATIC
      agrid->u[i][j][4] = p1/GAMM1
#ifdef MHD
	+ 0.5*(bx1*bx1 + by1*by1 + bz1*bz1) 
#endif
	+ 0.5*d1*(u1*u1 + v1*v1 + w1*w1);
#endif
    }
  }

  for(i=is1-2; i<im1; i++){
    for(j=im2; j<=ie2+2; j++){
      agrid->u[i][j][0] = d2;
      agrid->u[i][j][1] = u2*d2;
      agrid->u[i][j][2] = v2*d2;
      agrid->u[i][j][3] = w2*d2;
#ifdef MHD
      agrid->bx[i][j] = bx2;
      agrid->by[i][j] = by2;
      agrid->u[i][j][NVAR-3] = bx2;
      agrid->u[i][j][NVAR-2] = by2;
      agrid->u[i][j][NVAR-1] = bz2;
#endif
#ifdef ADIABATIC
      agrid->u[i][j][4] = p2/GAMM1 
#ifdef MHD
	+ 0.5*(bx2*bx2 + by2*by2 + bz2*bz2)
#endif
	+ 0.5*d2*(u2*u2 + v2*v2 + w2*w2);
#endif
    }
  }

  for(i=is1-2; i<im1; i++){
    for(j=is2-2; j<im2; j++){
      agrid->u[i][j][0] = d3;
      agrid->u[i][j][1] = u3*d3;
      agrid->u[i][j][2] = v3*d3;
      agrid->u[i][j][3] = w3*d3;
#ifdef MHD
      agrid->bx[i][j] = bx3;
      agrid->by[i][j] = by3;
      agrid->u[i][j][NVAR-3] = bx3;
      agrid->u[i][j][NVAR-2] = by3;
      agrid->u[i][j][NVAR-1] = bz3;
#endif
#ifdef ADIABATIC
      agrid->u[i][j][4] = p3/GAMM1 
#ifdef MHD
	+ 0.5*(bx3*bx3 + by3*by3 + bz3*bz3)
#endif
	+ 0.5*d3*(u3*u3 + v3*v3 + w3*w3);
#endif
    }
  }

  for(i=im1; i<=ie1+2; i++){
    for(j=is2-2; j<im2; j++){
      agrid->u[i][j][0] = d4;
      agrid->u[i][j][1] = u4*d4;
      agrid->u[i][j][2] = v4*d4;
      agrid->u[i][j][3] = w4*d4;
#ifdef MHD
      agrid->bx[i][j] = bx4;
      agrid->by[i][j] = by4;
      agrid->u[i][j][NVAR-3] = bx4;
      agrid->u[i][j][NVAR-2] = by4;
      agrid->u[i][j][NVAR-1] = bz4;
#endif
#ifdef ADIABATIC
      agrid->u[i][j][4] = p4/GAMM1 
#ifdef MHD
	+ 0.5*(bx4*bx4 + by4*by4 + bz4*bz4)
#endif
	+ 0.5*d4*(u4*u4 + v4*v4 + w4*w4);
#endif
    }
  }

  return;
}

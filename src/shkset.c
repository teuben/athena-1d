/*============================================================================*/
/*///////////////////////////// Function SHKSET \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
/*                                                                            */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

void shkset(FILE *p_input_file, struct grid_block *agrid)
{
/* Problem generator for 1-D Riemann problems                                 */
/*                                                                            */
/*============================================================================*/
  char buf[120];
  int i=0;
  int is,ie;
  Real dl,pl,ul,vl,wl,bxl,byl,bzl,dr,pr,ur,vr,wr,bxr,byr,bzr;
/*============================================================================*/

/* Read left state from 'athinput' */
/* Read dl,pl,ul,vl,wl,bxl,byl,bzl */

   if (fscanf(p_input_file,"%lg %[^\n]\n", &dl, buf) != 2) {
      printf("Error reading dl\n"); exit(EXIT_FAILURE);
   }
#ifdef ADIABATIC
   if (fscanf(p_input_file,"%lg %[^\n]\n", &pl, buf) != 2) {
      printf("Error reading pl\n"); exit(EXIT_FAILURE);
   }
#endif
   if (fscanf(p_input_file,"%lg %[^\n]\n", &ul, buf) != 2) {
      printf("Error reading ul\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &vl, buf) != 2) {
      printf("Error reading vl\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &wl, buf) != 2) {
      printf("Error reading wl\n"); exit(EXIT_FAILURE);
   }
#ifdef MHD
   if (fscanf(p_input_file,"%lg %[^\n]\n", &bxl, buf) != 2) {
      printf("Error reading bxl\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &byl, buf) != 2) {
      printf("Error reading byl\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &bzl, buf) != 2) {
      printf("Error reading bzl\n"); exit(EXIT_FAILURE);
   }
#endif

/* Read right state from 'athinput' */
/* Read dr,pr,ur,vr,wr,bxr,byr,bzr */

   if (fscanf(p_input_file,"%lg %[^\n]\n", &dr, buf) != 2) {
      printf("Error reading dr\n"); exit(EXIT_FAILURE);
   }
#ifdef ADIABATIC
   if (fscanf(p_input_file,"%lg %[^\n]\n", &pr, buf) != 2) {
      printf("Error reading pr\n"); exit(EXIT_FAILURE);
   }
#endif
   if (fscanf(p_input_file,"%lg %[^\n]\n", &ur, buf) != 2) {
      printf("Error reading ur\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &vr, buf) != 2) {
      printf("Error reading vr\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &wr, buf) != 2) {
      printf("Error reading wr\n"); exit(EXIT_FAILURE);
   }
#ifdef MHD
   if (fscanf(p_input_file,"%lg %[^\n]\n", &bxr, buf) != 2) {
      printf("Error reading bxr\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &byr, buf) != 2) {
      printf("Error reading byr\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &bzr, buf) != 2) {
      printf("Error reading bzr\n"); exit(EXIT_FAILURE);
   }
#endif

/* setup dependent variables for shocktube in X-direction */

   is = agrid->is; ie = agrid->ie;

   for (i=is-2; i<=is+((ie-is+1)/2)-1; i++) {
      agrid->u[i][0] = dl;
      agrid->u[i][1] = ul*dl;
      agrid->u[i][2] = vl*dl;
      agrid->u[i][3] = wl*dl;
#ifdef MHD
      agrid->bx[i] = bxl;
      agrid->u[i][NVAR-2] = byl;
      agrid->u[i][NVAR-1] = bzl;
#endif
#ifdef ADIABATIC
      agrid->u[i][4] = pl/GAMM1 
#ifdef MHD
        + 0.5*(bxl*bxl + byl*byl + bzl*bzl)
#endif
        + 0.5*dl*(ul*ul + vl*vl + wl*wl);
#endif
   }
   for (i=is+((ie-is+1)/2); i<=ie+2; i++) {
      agrid->u[i][0] = dr;
      agrid->u[i][1] = ur*dr;
      agrid->u[i][2] = vr*dr;
      agrid->u[i][3] = wr*dr;
#ifdef MHD
      agrid->bx[i] = bxr;
      agrid->u[i][NVAR-2] = byr;
      agrid->u[i][NVAR-1] = bzr;
#endif
#ifdef ADIABATIC
      agrid->u[i][4] = pr/GAMM1
#ifdef MHD
        + 0.5*(bxr*bxr + byr*byr + bzr*bzr) 
#endif
        + 0.5*dr*(ur*ur + vr*vr + wr*wr);
#endif
      }
}

/*============================================================================*/
/*///////////////////////// Function INIT_GRID_BLOCK \\\\\\\\\\\\\\\\\\\\\\\\\*/
/*                                                                            */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

void PROBLEM (FILE *p_input_file, struct grid_block *agrid);


#ifdef ONE_D
/* Wrap the function in an ifdef so that we can compile and link the
   1D and 2D init_grid_block files and only the appropriate one will
   contribute code.  Ultimately this ought to be handled at the
   configure step. */


void init_grid_block(FILE *p_input_file, FILE *p_output_file,
   struct grid_block *agrid)
{
/* Initializes data in grid_block structure agrid for new runs.  Calls        */
/* init_grid, PROBLEM, and init_dt                                            */
/*                                                                            */
/*============================================================================*/
  int i=0,nxzones;
  Real xmin,xmax;
  Real diib,sxiib,syiib,sziib,eiib;
  Real doib,sxoib,syoib,szoib,eoib;
#ifdef MHD
  Real bxiib,byiib,bziib;
  Real bxoib,byoib,bzoib;
#endif /* MHD */
  char buf[120];
/*============================================================================*/
/* Initialize grid */
/* Read lines 9-17 from 'athinput' */
/* Read nxzones,xmin,xmax */
/* Then initialize dx */

   if (fscanf(p_input_file,"%i %[^\n]\n", &nxzones, buf) != 2) {
      printf("Error reading nxzones in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
   } fprintf(p_output_file,"   NXZONES  = %i\n",nxzones);

   if (fscanf(p_input_file,"%lg %[^\n]\n", &xmin, buf) != 2) {
      printf("Error reading xmin in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
   } fprintf(p_output_file,"   XMIN     = %e\n",xmin);

   if (fscanf(p_input_file,"%lg %[^\n]\n", &xmax, buf) != 2) {
      printf("Error reading xmax in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
   } fprintf(p_output_file,"   XMAX     = %e\n",xmax);

   agrid->dx = (xmax-xmin)/(double)nxzones;

/* Test that number of zones in each dimension is OK */
/* Initialize is,ie */

   if (nxzones > NX1-8) {
      printf("Too many zones in X-direction, nxzones > NX1-8\n");
      exit(EXIT_FAILURE);
   } agrid->is = 4; agrid->ie=nxzones+3;

/* Read boundary condition flags */

   if (fscanf(p_input_file,"%i %[^\n]\n", &agrid->niib, buf) != 2) {
      printf("Error reading niib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
   } fprintf(p_output_file,"   NIIB     = %i\n",agrid->niib);

   if (fscanf(p_input_file,"%i %[^\n]\n", &agrid->noib, buf) != 2) {
      printf("Error reading noib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
   } fprintf(p_output_file,"   NOIB     = %i\n",agrid->noib);

/* If inner-I boundary is flow-in, read constant states */

   if (agrid->niib == 3) {
     if (fscanf(p_input_file,"%lg %[^\n]\n",&diib,buf) != 2) {
        printf("Error reading diib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   DIIB     = %e\n",diib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&sxiib,buf) != 2) {
        printf("Error reading sxiib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   SXIIB    = %e\n",sxiib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&syiib,buf) != 2) {
        printf("Error reading syiib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   SYIIB    = %e\n",syiib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&sziib,buf) != 2) {
        printf("Error reading sziib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   SZIIB    = %e\n",sziib);
#ifdef ADIABATIC
     if (fscanf(p_input_file,"%lg %[^\n]\n",&eiib,buf) != 2) {
        printf("Error reading eiib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   EIIB     = %e\n",eiib);
#endif
#ifdef MHD
     if (fscanf(p_input_file,"%lg %[^\n]\n",&bxiib,buf) != 2) {
        printf("Error reading bxiib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   BXIIB    = %e\n",bxiib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&byiib,buf) != 2) {
        printf("Error reading byiib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   BYIIB    = %e\n",byiib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&bziib,buf) != 2) {
        printf("Error reading bziib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   BZIIB    = %e\n",bziib);
#endif
     for (i=0; i<4;   i++) {
       agrid->boundary_values.uiib[i][0] = diib;
       agrid->boundary_values.uiib[i][1] = sxiib;
       agrid->boundary_values.uiib[i][2] = syiib;
       agrid->boundary_values.uiib[i][3] = sziib;
#ifdef ADIABATIC
       agrid->boundary_values.uiib[i][4] = eiib;
#endif
#ifdef MHD
       agrid->boundary_values.bxiib[i] = bxiib;
       agrid->boundary_values.uiib[i][NVAR-3] = bxiib;
       agrid->boundary_values.uiib[i][NVAR-2] = byiib;
       agrid->boundary_values.uiib[i][NVAR-1] = bziib;
#endif
     }
   }

/* If outer-I boundary is flow-in, read constant states */

   if (agrid->noib == 3) {
     if (fscanf(p_input_file,"%lg %[^\n]\n",&doib,buf) != 2) {
        printf("Error reading doib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   DOIB     = %e\n",doib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&sxoib,buf) != 2) {
        printf("Error reading sxoib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   SXOIB    = %e\n",sxoib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&syoib,buf) != 2) {
        printf("Error reading syoib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   SYOIB    = %e\n",syoib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&szoib,buf) != 2) {
        printf("Error reading szoib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   SZOIB    = %e\n",szoib);
#ifdef ADIABATIC
     if (fscanf(p_input_file,"%lg %[^\n]\n",&eoib,buf) != 2) {
        printf("Error reading eoib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   EOIB     = %e\n",eoib);
#endif
#ifdef MHD
     if (fscanf(p_input_file,"%lg %[^\n]\n",&bxoib,buf) != 2) {
        printf("Error reading bxoib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   BXOIB    = %e\n",bxoib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&byoib,buf) != 2) {
        printf("Error reading byoib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   BYOIB    = %e\n",byoib);

     if (fscanf(p_input_file,"%lg %[^\n]\n",&bzoib,buf) != 2) {
        printf("Error reading bzoib in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
     } fprintf(p_output_file,"   BZOIB    = %e\n",bzoib);
#endif
     for (i=0; i<4;   i++) {
       agrid->boundary_values.uoib[i][0] = doib;
       agrid->boundary_values.uoib[i][1] = sxoib;
       agrid->boundary_values.uoib[i][2] = syoib;
       agrid->boundary_values.uoib[i][3] = szoib;
#ifdef ADIABATIC
       agrid->boundary_values.uoib[i][4] = eoib;
#endif
#ifdef MHD
       agrid->boundary_values.bxoib[i] = bxoib;
       agrid->boundary_values.uoib[i][NVAR-3] = bxoib;
       agrid->boundary_values.uoib[i][NVAR-2] = byoib;
       agrid->boundary_values.uoib[i][NVAR-1] = bzoib;
#endif
     }
   }

/* Initialize data (in user-supplied routine PROBLEM), initialize dt */

   PROBLEM (p_input_file, agrid);
   agrid->time = 0.0;
   agrid->dt = init_dt(agrid);
}

#endif /* ONE_D */

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


#ifdef TWO_D
/* Wrap the function in an ifdef so that we can compile and link the
   1D and 2D init_grid_block files and only the appropriate one will
   contribute code.  Ultimately this ought to be handled at the
   configure step. */


void init_grid_block(FILE *p_input_file,FILE *p_output_file,Grid_Block *agrid){
/* Initializes data in grid_block structure agrid for new runs.  Calls        */
/* init_grid, PROBLEM, and init_dt                                            */
/*                                                                            */
/*============================================================================*/
  int i,j,nxzones,nyzones;
  Real xmin,xmax,ymin,ymax;
  Real diib1,sxiib1,syiib1,sziib1,eiib1,bxiib1,byiib1,bziib1;
  Real doib1,sxoib1,syoib1,szoib1,eoib1,bxoib1,byoib1,bzoib1;
  Real diib2,sxiib2,syiib2,sziib2,eiib2,bxiib2,byiib2,bziib2;
  Real doib2,sxoib2,syoib2,szoib2,eoib2,bxoib2,byoib2,bzoib2;
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
  } agrid->is1 = 4; agrid->ie1 = nxzones+3;

/* Read boundary condition flags */

  if (fscanf(p_input_file,"%i %[^\n]\n", &agrid->niib1, buf) != 2) {
    printf("Error reading niib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   NIIB1     = %i\n",agrid->niib1);

  if (fscanf(p_input_file,"%i %[^\n]\n", &agrid->noib1, buf) != 2) {
    printf("Error reading noib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   NOIB1     = %i\n",agrid->noib1);

/*============================================================================*/

  if (fscanf(p_input_file,"%i %[^\n]\n", &nyzones, buf) != 2) {
    printf("Error reading nyzones in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   NYZONES  = %i\n",nyzones);

  if (fscanf(p_input_file,"%lg %[^\n]\n", &ymin, buf) != 2) {
    printf("Error reading ymin in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   YMIN     = %e\n",ymin);

  if (fscanf(p_input_file,"%lg %[^\n]\n", &ymax, buf) != 2) {
    printf("Error reading ymax in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   YMAX     = %e\n",ymax);

  agrid->dy = (ymax-ymin)/(double)nyzones;

/* Test that number of zones in each dimension is OK */
/* Initialize is,ie */

  if (nyzones > NX2-8) {
    printf("Too many zones in X-direction, nyzones > NX2-8\n");
    exit(EXIT_FAILURE);
  } agrid->is2 = 4; agrid->ie2 = nyzones+3;

/* Read boundary condition flags */

  if (fscanf(p_input_file,"%i %[^\n]\n", &agrid->niib2, buf) != 2) {
    printf("Error reading niib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   NIIB2     = %i\n",agrid->niib2);

  if (fscanf(p_input_file,"%i %[^\n]\n", &agrid->noib2, buf) != 2) {
    printf("Error reading noib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   NOIB2     = %i\n",agrid->noib2);


/*============================================================================*/


/* If inner-1 boundary is flow-in, read constant states */

  if (agrid->niib1 == 3) {
    if (fscanf(p_input_file,"%lg %[^\n]\n",&diib1,buf) != 2) {
      printf("Error reading diib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   DIIB1     = %e\n",diib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&sxiib1,buf) != 2) {
      printf("Error reading sxiib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SXIIB1    = %e\n",sxiib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&syiib1,buf) != 2) {
      printf("Error reading syiib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SYIIB1    = %e\n",syiib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&sziib1,buf) != 2) {
      printf("Error reading sziib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SZIIB1    = %e\n",sziib1);
#ifdef ADIABATIC
    if (fscanf(p_input_file,"%lg %[^\n]\n",&eiib1,buf) != 2) {
      printf("Error reading eiib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   EIIB1     = %e\n",eiib1);
#endif
#ifdef MHD
    if (fscanf(p_input_file,"%lg %[^\n]\n",&bxiib1,buf) != 2) {
      printf("Error reading bxiib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BXIIB1    = %e\n",bxiib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&byiib1,buf) != 2) {
      printf("Error reading byiib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BYIIB1    = %e\n",byiib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&bziib1,buf) != 2) {
      printf("Error reading bziib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BZIIB1    = %e\n",bziib1);
#endif
    for(i=0; i<4; i++){
      for(j=agrid->is2-4; j<=agrid->ie2+4; j++){
	agrid->boundary_values.uiib1[i][j][0] = diib1;
	agrid->boundary_values.uiib1[i][j][1] = sxiib1;
	agrid->boundary_values.uiib1[i][j][2] = syiib1;
	agrid->boundary_values.uiib1[i][j][3] = sziib1;
#ifdef ADIABATIC
	agrid->boundary_values.uiib1[i][j][4] = eiib1;
#endif
#ifdef MHD
	agrid->boundary_values.bxiib1[i][j] = bxiib1;
	agrid->boundary_values.byiib1[i][j] = byiib1;
	agrid->boundary_values.uiib1[i][j][NVAR-3] = bxiib1;
	agrid->boundary_values.uiib1[i][j][NVAR-2] = byiib1;
	agrid->boundary_values.uiib1[i][j][NVAR-1] = bziib1;
#endif
      }
    }
  }

/* If outer-1 boundary is flow-in, read constant states */

  if (agrid->noib1 == 3) {
    if (fscanf(p_input_file,"%lg %[^\n]\n",&doib1,buf) != 2) {
      printf("Error reading doib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   DOIB1     = %e\n",doib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&sxoib1,buf) != 2) {
      printf("Error reading sxoib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SXOIB1    = %e\n",sxoib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&syoib1,buf) != 2) {
      printf("Error reading syoib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SYOIB1    = %e\n",syoib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&szoib1,buf) != 2) {
      printf("Error reading szoib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SZOIB1    = %e\n",szoib1);
#ifdef ADIABATIC
    if (fscanf(p_input_file,"%lg %[^\n]\n",&eoib1,buf) != 2) {
      printf("Error reading eoib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   EOIB1     = %e\n",eoib1);
#endif
#ifdef MHD
    if (fscanf(p_input_file,"%lg %[^\n]\n",&bxoib1,buf) != 2) {
      printf("Error reading bxoib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BXOIB1    = %e\n",bxoib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&byoib1,buf) != 2) {
      printf("Error reading byoib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BYOIB1    = %e\n",byoib1);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&bzoib1,buf) != 2) {
      printf("Error reading bzoib1 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BZOIB1    = %e\n",bzoib1);
#endif
    for(i=0; i<4; i++){
      for(j=agrid->is2-4; j<=agrid->ie2+4; j++){
	agrid->boundary_values.uoib1[i][j][0] = doib1;
	agrid->boundary_values.uoib1[i][j][1] = sxoib1;
	agrid->boundary_values.uoib1[i][j][2] = syoib1;
	agrid->boundary_values.uoib1[i][j][3] = szoib1;
#ifdef ADIABATIC
	agrid->boundary_values.uoib1[i][j][4] = eoib1;
#endif
#ifdef MHD
	agrid->boundary_values.bxoib1[i][j] = bxoib1;
	agrid->boundary_values.byoib1[i][j] = byoib1;
	agrid->boundary_values.uoib1[i][j][NVAR-3] = bxoib1;
	agrid->boundary_values.uoib1[i][j][NVAR-2] = byoib1;
	agrid->boundary_values.uoib1[i][j][NVAR-1] = bzoib1;
#endif
      }
    }
  }

/*============================================================================*/

/* If inner-2 boundary is flow-in, read constant states */

  if (agrid->niib2 == 3) {
    if (fscanf(p_input_file,"%lg %[^\n]\n",&diib2,buf) != 2) {
      printf("Error reading diib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   DIIB2     = %e\n",diib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&sxiib2,buf) != 2) {
      printf("Error reading sxiib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SXIIB2    = %e\n",sxiib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&syiib2,buf) != 2) {
      printf("Error reading syiib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SYIIB2    = %e\n",syiib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&sziib2,buf) != 2) {
      printf("Error reading sziib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SZIIB2    = %e\n",sziib2);
#ifdef ADIABATIC
    if (fscanf(p_input_file,"%lg %[^\n]\n",&eiib2,buf) != 2) {
      printf("Error reading eiib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   EIIB2     = %e\n",eiib2);
#endif
#ifdef MHD
    if (fscanf(p_input_file,"%lg %[^\n]\n",&bxiib2,buf) != 2) {
      printf("Error reading bxiib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BXIIB2    = %e\n",bxiib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&byiib2,buf) != 2) {
      printf("Error reading byiib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BYIIB2    = %e\n",byiib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&bziib2,buf) != 2) {
      printf("Error reading bziib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BZIIB2    = %e\n",bziib2);
#endif
    for(i=agrid->is1-4; i<=agrid->ie1+4; i++){
      for(j=0; j<4; j++){
	agrid->boundary_values.uiib2[i][j][0] = diib2;
	agrid->boundary_values.uiib2[i][j][1] = sxiib2;
	agrid->boundary_values.uiib2[i][j][2] = syiib2;
	agrid->boundary_values.uiib2[i][j][3] = sziib2;
#ifdef ADIABATIC
	agrid->boundary_values.uiib2[i][j][4] = eiib2;
#endif
#ifdef MHD
	agrid->boundary_values.bxiib2[i][j] = bxiib2;
	agrid->boundary_values.byiib2[i][j] = byiib2;
	agrid->boundary_values.uiib2[i][j][NVAR-3] = bxiib2;
	agrid->boundary_values.uiib2[i][j][NVAR-2] = byiib2;
	agrid->boundary_values.uiib2[i][j][NVAR-1] = bziib2;
#endif
      }
    }
  }

/* If outer-2 boundary is flow-in, read constant states */

  if (agrid->noib2 == 3) {
    if (fscanf(p_input_file,"%lg %[^\n]\n",&doib2,buf) != 2) {
      printf("Error reading doib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   DOIB2     = %e\n",doib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&sxoib2,buf) != 2) {
      printf("Error reading sxoib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SXOIB2    = %e\n",sxoib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&syoib2,buf) != 2) {
      printf("Error reading syoib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SYOIB2    = %e\n",syoib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&szoib2,buf) != 2) {
      printf("Error reading szoib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   SZOIB2    = %e\n",szoib2);
#ifdef ADIABATIC
    if (fscanf(p_input_file,"%lg %[^\n]\n",&eoib2,buf) != 2) {
      printf("Error reading eoib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   EOIB2     = %e\n",eoib2);
#endif
#ifdef MHD
    if (fscanf(p_input_file,"%lg %[^\n]\n",&bxoib2,buf) != 2) {
      printf("Error reading bxoib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BXOIB2    = %e\n",bxoib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&byoib2,buf) != 2) {
      printf("Error reading byoib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BYOIB2    = %e\n",byoib2);

    if (fscanf(p_input_file,"%lg %[^\n]\n",&bzoib2,buf) != 2) {
      printf("Error reading bzoib2 in INIT_GRID_BLOCK\n"); exit(EXIT_FAILURE);
    } fprintf(p_output_file,"   BZOIB2    = %e\n",bzoib2);
#endif
    for(i=agrid->is1-4; i<=agrid->ie1+4; i++){
      for(j=0; j<4; j++){
	agrid->boundary_values.uoib2[i][j][0] = doib2;
	agrid->boundary_values.uoib2[i][j][1] = sxoib2;
	agrid->boundary_values.uoib2[i][j][2] = syoib2;
	agrid->boundary_values.uoib2[i][j][3] = szoib2;
#ifdef ADIABATIC
	agrid->boundary_values.uoib2[i][j][4] = eoib2;
#endif
#ifdef MHD
	agrid->boundary_values.bxoib2[i][j] = bxoib2;
	agrid->boundary_values.byoib2[i][j] = byoib2;
	agrid->boundary_values.uoib2[i][j][NVAR-3] = bxoib2;
	agrid->boundary_values.uoib2[i][j][NVAR-2] = byoib2;
	agrid->boundary_values.uoib2[i][j][NVAR-1] = bzoib2;
#endif
      }
    }
  }

/* Initialize data (in user-supplied routine PROBLEM), initialize dt */

  PROBLEM (p_input_file, agrid);
  agrid->time = 0.0;
  agrid->dt = init_dt(agrid);
}

#endif /* TWO_D */

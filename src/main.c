/*============================================================================*/
/*/////////////////////////// ATHENA Main Program \\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
/*                                                                            */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

int main(void)
{
/* Temporary driver program to test Roe's method.                             */
/*                                                                            */
/*============================================================================*/
  struct grid_block grid_level0;
  struct bval_array bval_level0;
  int ires,nlim;
  Real t_res,dt_res,tlim,cpu_time,zcs;
  char buf[120],id[3],res_file[9];
  FILE *p_input_file,*p_output_file;
  clock_t time0,time1;
  time_t start;
/*============================================================================*/
/* Open and read job control parameters from lines 1-6 of 'athinput' */

  if((p_input_file = fopen("athinput", "r")) == NULL){
    fprintf(stderr,"[main]: Unable to open \"athinput\" file\n");
    exit(EXIT_FAILURE);
  }
  if((p_output_file = fopen("athoutput", "wb")) == NULL){
    fprintf(stderr,"[main]: Unable to open \"athoutput\" file\n");
    fclose(p_input_file);
    exit(EXIT_FAILURE);
  }
  if(time(&start)>0){ /* current calendar time (UTC) is available */
    fprintf(p_output_file,"OUTPUT FILE for ATHENA created on %s\n",ctime(&start));
  }
  else{
    fprintf(p_output_file,"OUTPUT FILE for ATHENA created on XXX at XXX\n\n");
  }

  if (fscanf(p_input_file,"%s %[^\n]\n", id, buf) != 2) {
    printf("Error reading id in MAIN\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   ID       = %s\n",id);

  if (fscanf(p_input_file,"%i %[^\n]\n", &ires, buf) != 2) {
    printf("Error reading ires in MAIN\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   IRES     = %i\n",ires);

  if (fscanf(p_input_file,"%s %[^\n]\n", res_file, buf) != 2) {
    printf("Error reading res_file in MAIN\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   RES_FILE = %s\n",res_file);

  if (fscanf(p_input_file,"%lg %[^\n]\n", &dt_res, buf) != 2) {
    printf("Error reading dt_res in MAIN\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   DT_RES   = %e\n",dt_res);

  if (fscanf(p_input_file,"%i %[^\n]\n", &nlim, buf) != 2) {
    printf("Error reading nlim in MAIN\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   NLIM     = %i\n",nlim);

  if (fscanf(p_input_file,"%lg %[^\n]\n", &tlim, buf) != 2) {
    printf("Error reading tlim in MAIN\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   TLIM     = %e\n",tlim);

/* Initialize grid_block level0. */

  if (ires == 1) {
/*      mget(res_file, &grid_level0); */
    t_res = grid_level0.time;
  }

  if (fscanf(p_input_file,"%lg %[^\n]\n", &grid_level0.dt_hst, buf) != 2) {
    printf("Error reading dt_hst in MAIN\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   DT_HST   = %e\n",grid_level0.dt_hst);

  if (fscanf(p_input_file,"%lg %[^\n]\n", &grid_level0.dt_hdf, buf) != 2) {
    printf("Error reading dt_hdf in MAIN\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   DT_HDF   = %e\n",grid_level0.dt_hdf);

  if (fscanf(p_input_file,"%lg %[^\n]\n", &grid_level0.dt_bin, buf) != 2) {
    printf("Error reading dt_bin in MAIN\n"); exit(EXIT_FAILURE);
  } fprintf(p_output_file,"   DT_BIN   = %e\n",grid_level0.dt_bin);
  printf("dt_bin = %e\n",grid_level0.dt_bin);

  if (ires == 0) {
    init_grid_block(p_input_file, p_output_file, &grid_level0);
    t_res = 0.0;
    grid_level0.t_hst = 0.0;
    grid_level0.t_bin = 0.0;
    grid_level0.t_hdf = 0.0;
    grid_level0.ncycles = 0;
    sprintf(grid_level0.hst_file, "hst000%2s", id);
    sprintf(grid_level0.hdf_file, "hdf000%2s", id);
    sprintf(grid_level0.bin_file, "bin000%2s", id);
  }

/* Setup complete, dump initial conditions */

  if (ires == 0) {

    set_bval_arrays(&grid_level0, &bval_level0); 
    set_ghost_zones(&grid_level0, &bval_level0); 

    data_output(&grid_level0, 1, 1, 1);
  }
  fprintf(p_output_file,"\nSetup complete, entering main loop...\n\n");
  printf("\nSetup complete, entering main loop...\n\n");
  fclose(p_input_file);

/* Integrate until stopping criteria reached */

  time0 = clock();
  while (grid_level0.ncycles < nlim AND grid_level0.time < tlim) {

    set_bval_arrays(&grid_level0, &bval_level0); 
    set_ghost_zones(&grid_level0, &bval_level0); 

    if ((tlim-grid_level0.time) < grid_level0.dt) 
      {grid_level0.dt=(tlim-grid_level0.time);}
    integrate_1step(&grid_level0);

    data_output(&grid_level0, 0, 0, 0);

    printf("cycle=%i time=%e dt=%e\n",grid_level0.ncycles,grid_level0.time,
	   grid_level0.dt);

  }

/*********************/
  linear_wave_error(&grid_level0);
/*********************/

/* Finish up */

  time1 = clock();
  cpu_time = (double)(time1-time0)/(double)CLOCKS_PER_SEC;
  zcs = (grid_level0.ncycles)*(grid_level0.ie-grid_level0.is+1)/cpu_time;

  if (grid_level0.ncycles == nlim) {
    printf("\nterminating on cycle limit\n");
  } else {
    printf("\nterminating on time limit\n");
  }
  printf("  tlim= %e   nlim= %i\n",tlim,nlim);
  printf("  time= %e  cycle= %i\n",grid_level0.time,grid_level0.ncycles);
  printf("\nzone-cycles/cpu-second = %e\n\n",zcs);
  set_bval_arrays(&grid_level0, &bval_level0); 
  set_ghost_zones(&grid_level0, &bval_level0); 
  data_output(&grid_level0, 1, 1, 1);
  fclose(p_output_file);
}

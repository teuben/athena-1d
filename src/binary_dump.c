/*==============================================================================
 *                          Function BINARY_DUMP
 * Writes an unformatted dump of the field variables to files named BINXXXID
 * where XXX is a 3-digit integer, and ID is a two-letter problem tag.
 * Normally only one of BINARY or HDF file dumps are used to save data for
 * analysis.
 *
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

void binary_dump(struct grid_block *agrid)
{
  FILE *p_binfile;
  int ndata[2],nzones,is,ie,i,n;
  float eos[2],*data,*xgrid;
  is = agrid->is;
  ie = agrid->ie;
  p_binfile = fopen(agrid->bin_file,"wb");

  /* Write number of zones and variables */

  ndata[0] = ie-is+1;
  ndata[1] = NVAR;
#ifdef MHD
  ndata[1] = NVAR + 1;
#endif
  fwrite(ndata,sizeof(int),2,p_binfile);

  /* Write (gamma-1) and isothermal sound speed */

  eos[0] = (float)GAMM1 ;
  eos[1] = (float)ISOTHERMAL_C ;
  fwrite(eos,sizeof(float),2,p_binfile);

  /* Write grid */

  nzones = ndata[0];
  xgrid = (float *) malloc(ndata[0]*sizeof(float));
  xgrid[0] = 0.5*(float)agrid->dx;
  for (i=1;i<ndata[0];i++) xgrid[i] = xgrid[i-1] + (float)agrid->dx;
  fwrite(xgrid,sizeof(float),ndata[0],p_binfile);

  /* Write data */

  data = (float *) malloc(ndata[0]*sizeof(float));
  for (n=0;n<NVAR; n++) {
    for (i=0;i<ndata[0]; i++) {data[i] = (float)agrid->u[i+is][n];}
    fwrite(data,sizeof(float),ndata[0],p_binfile);
  }

  for (i=0;i<ndata[0]; i++) {data[i] = (float)agrid->bx[i+is];}
  fwrite(data,sizeof(float),ndata[0],p_binfile);

  /* Close dump file, increment filename */

  fclose(p_binfile);
  add1_2name(agrid->bin_file);
}

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
#include <string.h>
#include <endian.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

#ifdef WRITE_DX_HEADER
static void write_dx_header(const struct grid_block *agrid){
  int n; /* Dummy loop var */
  int offset; /* Byte offset */
  char FileName[20]; /* Assume 20 char. is sufficient */
  FILE *pfile; /* Pointer to file. */
  int nzones;

  strcpy(FileName,agrid->bin_file);/* Copy the Data Dump File Name */
  strcat(FileName,".dx"); /* Append ".dx" for the header name */
  nzones = agrid->ie - agrid->is + 1;

  /* Open the output file with the name "FileName" just constructed. */
  if ((pfile = fopen(FileName,"w")) == NULL) {
    fprintf(stderr,"[write_dx_header]: File Open Error Occured\n");
    return; /* Error */
  }

  fprintf(pfile,"# The default data type is ");
#if __BYTE_ORDER == __LITTLE_ENDIAN
  fprintf(pfile,"least-significant-byte first and binary\n");
  fprintf(pfile,"data mode lsb binary\n\n");
#elif __BYTE_ORDER == __BIG_ENDIAN
  fprintf(pfile,"most-significant-byte first and binary\n");
  fprintf(pfile,"data mode msb binary\n\n");
#else
#error : Machine is neither BIG ENDIAN or LITTLE ENDIAN
#endif /* __BYTE_ORDER */

  fprintf(pfile,"# object 1 is the regular positions\n");
  fprintf(pfile,"object 1 class gridpositions counts %d\n",nzones);
  fprintf(pfile,"origin %g\n",0.5*agrid->dx);
  fprintf(pfile,"delta  %g\n\n",agrid->dx);

  fprintf(pfile,"# objects 2 on are the data, which ");
  fprintf(pfile,"are in a one-to-one correspondence\n");
  fprintf(pfile,"# with the positions (\"dep\" on positions).\n\n");

  /* Initialize the byte offset to skip over the first 4 floats */
  offset = 4*sizeof(float);

  offset += nzones*sizeof(float);
  fprintf(pfile,"# mass density\n");
  fprintf(pfile,"object 2 class array type float rank 0 items %d\n",nzones);
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");

  offset += nzones*sizeof(float);
  fprintf(pfile,"# 1-momentum\n");
  fprintf(pfile,"object 3 class array type float rank 0 items %d\n",nzones);
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");

  offset += nzones*sizeof(float);
  fprintf(pfile,"# 2-momentum\n");
  fprintf(pfile,"object 4 class array type float rank 0 items %d\n",nzones);
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");

  offset += nzones*sizeof(float);
  fprintf(pfile,"# 3-momentum\n");
  fprintf(pfile,"object 5 class array type float rank 0 items %d\n",nzones);
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");

  offset += nzones*sizeof(float);
  fprintf(pfile,"# total energy density\n");
  fprintf(pfile,"object 6 class array type float rank 0 items %d\n",nzones);
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");

#ifdef MHD
  offset += nzones*sizeof(float);
  fprintf(pfile,"# 2-field\n");
  fprintf(pfile,"object 7 class array type float rank 0 items %d\n",nzones);
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");

  offset += nzones*sizeof(float);
  fprintf(pfile,"# 3-field\n");
  fprintf(pfile,"object 8 class array type float rank 0 items %d\n",nzones);
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");

  offset += nzones*sizeof(float);
  fprintf(pfile,"# 1-field\n");
  fprintf(pfile,"object 9 class array type float rank 0 items %d\n",nzones);
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");
#endif /* MHD */

  for(n=0;n<72;n++) fputc((int)'#',pfile);

  fprintf(pfile,"\n\n# Next we construct a field for each of the objects above.\n");

  fprintf(pfile,"object \"Mass Density\" class field\n");
  fprintf(pfile,"component \"positions\" value 1\n");
  fprintf(pfile,"component \"data\" value 2\n\n");

  fprintf(pfile,"object \"1-Momentum\" class field\n");
  fprintf(pfile,"component \"positions\" value 1\n");
  fprintf(pfile,"component \"data\" value 3\n\n");

  fprintf(pfile,"object \"2-Momentum\" class field\n");
  fprintf(pfile,"component \"positions\" value 1\n");
  fprintf(pfile,"component \"data\" value 4\n\n");

  fprintf(pfile,"object \"3-Momentum\" class field\n");
  fprintf(pfile,"component \"positions\" value 1\n");
  fprintf(pfile,"component \"data\" value 5\n\n");

  fprintf(pfile,"object \"Energy Density\" class field\n");
  fprintf(pfile,"component \"positions\" value 1\n");
  fprintf(pfile,"component \"data\" value 6\n\n");

#ifdef MHD
  fprintf(pfile,"object \"2-Field\" class field\n");
  fprintf(pfile,"component \"positions\" value 1\n");
  fprintf(pfile,"component \"data\" value 7\n\n");

  fprintf(pfile,"object \"3-Field\" class field\n");
  fprintf(pfile,"component \"positions\" value 1\n");
  fprintf(pfile,"component \"data\" value 8\n\n");

  fprintf(pfile,"object \"1-Field\" class field\n");
  fprintf(pfile,"component \"positions\" value 1\n");
  fprintf(pfile,"component \"data\" value 9\n\n");
#endif /* MHD */

  for(n=0;n<72;n++) fputc((int)'#',pfile);

  fprintf(pfile,"\n\n# Now we construct a string object");
  fprintf(pfile," which stores the names of the fields.\n");
  fprintf(pfile,"# This way we can select an individual field");
  fprintf(pfile," before importing the data.\n\n");

  fprintf(pfile,"object \"Variables\" class string\n");
  fprintf(pfile," \"Mass Density\"\n");
  fprintf(pfile," \"1-Momentum\"\n");
  fprintf(pfile," \"2-Momentum\"\n");
  fprintf(pfile," \"3-Momentum\"\n");
#ifdef MHD
  fprintf(pfile," \"1-Field\"\n");
  fprintf(pfile," \"2-Field\"\n");
  fprintf(pfile," \"3-Field\"\n");
#endif /* MHD */
  fprintf(pfile," \"Energy Density\"\n");

  fprintf(pfile,"\nend\n");

  fclose(pfile);
  return;
}
#endif /* WRITE_DX_HEADER */



void binary_dump(struct grid_block *agrid)
{
  FILE *p_binfile;
  int ndata[2],nzones,is,ie,i,n;
  float eos[2],*data,*xgrid;

  if((p_binfile = fopen(agrid->bin_file,"wb")) == NULL){
    fprintf(stderr,"[binary_dump]: Unable to open \"%s\" file\n",
	    agrid->bin_file);
    return;
  }
  is = agrid->is;
  ie = agrid->ie;
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

#ifdef MHD
  for (i=0;i<ndata[0]; i++) {data[i] = (float)agrid->bx[i+is];}
  fwrite(data,sizeof(float),ndata[0],p_binfile);
#endif

  /* Close dump file, increment filename */

  fclose(p_binfile);
#ifdef WRITE_DX_HEADER
  write_dx_header(agrid);
#endif /* WRITE_DX_HEADER */
  add1_2name(agrid->bin_file);
}

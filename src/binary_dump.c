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
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

#ifdef ONE_D
/* Wrap the functions in an ifdef so that we can compile and link the
   1D and 2D binary data dump files and only the appropriate one will
   contribute code.  Ultimately this ought to be handled at the
   configure step. */

#ifdef WRITE_DX_HEADER

#ifdef MHD
#define NUM_ARRAY (NVAR + 1)
#else
#define NUM_ARRAY (NVAR)
#endif

static int big_endian(){
  short int n = 1;
  char *ep = (char *)&n;

  return (*ep == 0); /* Returns 1 on a big endian machine */
}



static void write_dx_header(const struct grid_block *agrid){
  int n; /* Dummy loop var */
  int offset; /* Byte offset */
  char FileName[20]; /* Assume 20 char. is sufficient */
  FILE *pfile; /* Pointer to file. */
  int nzones;
  char array_name[][18] = 
    {"Mass Density",
     "1-Momentum",
     "2-Momentum",
     "3-Momentum"
#ifdef ADIABATIC
     ,"Energy Density"
#endif
#ifdef MHD
     ,"1-Field","2-Field","3-Field","1-Interface-Field"
#endif
    };

  strcpy(FileName,agrid->bin_file);/* Copy the Data Dump File Name */
  strcat(FileName,".dx"); /* Append ".dx" for the header name */
  nzones = agrid->ie - agrid->is + 1;

  /* Open the output file with the name "FileName" just constructed. */
  if ((pfile = fopen(FileName,"w")) == NULL) {
    fprintf(stderr,"[write_dx_header]: File Open Error Occured\n");
    return; /* Error */
  }

  fprintf(pfile,"# The default data type is ");
  if(big_endian()){
    fprintf(pfile,"most-significant-byte first and binary\n");
    fprintf(pfile,"data mode msb binary\n\n");
  }
  else{
    fprintf(pfile,"least-significant-byte first and binary\n");
    fprintf(pfile,"data mode lsb binary\n\n");
  }

  fprintf(pfile,"# object 1 is the cell center positions\n");
  fprintf(pfile,"object 1 class gridpositions counts %d\n",nzones);
  fprintf(pfile,"origin %g\n",0.5*agrid->dx);
  fprintf(pfile,"delta  %g\n\n",agrid->dx);

  fprintf(pfile,"# object 2 is the cell interface positions\n");
  fprintf(pfile,"object 2 class gridpositions counts %d\n",nzones+1);
  fprintf(pfile,"origin %g\n",0.0);
  fprintf(pfile,"delta  %g\n\n",agrid->dx);

  fprintf(pfile,"# objects 3 on are the data, which ");
  fprintf(pfile,"are in a one-to-one correspondence\n");
  fprintf(pfile,"# with the positions (\"dep\" on positions).\n\n");

  /* Initialize the byte offset to skip over the first 4 floats */
  offset = 2*sizeof(int) + 2*sizeof(float);

  for(n=0;n<NUM_ARRAY;n++){
    offset += nzones*sizeof(float);
    fprintf(pfile,"# %s\n",array_name[n]);
#ifdef MHD
    /* The last item is the x-interface field for which there is 1
       extra element */
    fprintf(pfile,"object %d class array type float rank 0 items %d\n",
	    n+3,(n == (NUM_ARRAY - 1) ? nzones + 1 : nzones));
#else
    fprintf(pfile,"object %d class array type float rank 0 items %d\n",
	    n+3,nzones);
#endif
    fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
    fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");
  }

  for(n=0;n<72;n++) fputc((int)'#',pfile);

  fprintf(pfile,"\n\n# Next we construct a field for each of the objects above.\n");

  for(n=0;n<NUM_ARRAY;n++){
    fprintf(pfile,"object \"%s\" class field\n",array_name[n]);
    fprintf(pfile,"component \"positions\" value %d\n",
	    (n == (NUM_ARRAY - 1) ? 2 : 1));
    fprintf(pfile,"component \"data\" value %d\n\n",n+3);
  }

  for(n=0;n<72;n++) fputc((int)'#',pfile);

  fprintf(pfile,"\n\n# Now we construct a string object");
  fprintf(pfile," which stores the names of the fields.\n");
  fprintf(pfile,"# This way we can select an individual field");
  fprintf(pfile," before importing the data.\n\n");

  fprintf(pfile,"object \"Variables\" class string\n");
  for(n=0;n<NUM_ARRAY;n++){
    fprintf(pfile," \"%s\"\n",array_name[n]);
  }

  fprintf(pfile,"\nend\n");

  fclose(pfile);
  return;
}
#undef NUM_ARRAY
#endif /* WRITE_DX_HEADER */



void binary_dump(struct grid_block *agrid)
{
  FILE *p_binfile;
  int ndata[2],is,ie,i,n;
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
  fwrite(ndata,sizeof(int),2,p_binfile);

  /* Write (gamma-1) and isothermal sound speed */

  eos[0] = (float)GAMM1 ;
  eos[1] = (float)ISOTHERMAL_C ;
  fwrite(eos,sizeof(float),2,p_binfile);

  /* Write grid */

  if((xgrid = (float *)malloc((1+ndata[0])*sizeof(float))) == NULL){
    fprintf(stderr,"[binary_dump]: Unable to malloc memory to write binary file %s\n",agrid->bin_file);
    fclose(p_binfile);
    add1_2name(agrid->bin_file);
    return;
  }

  xgrid[0] = 0.5*(float)agrid->dx;
  for (i=1;i<ndata[0];i++) xgrid[i] = xgrid[i-1] + (float)agrid->dx;
  fwrite(xgrid,sizeof(float),ndata[0],p_binfile);

  /* Write data */

  data = xgrid; /* Copy a pointer to the memory we just malloc'd */
  for (n=0;n<NVAR; n++) {
    for (i=0;i<ndata[0]; i++) {data[i] = (float)agrid->u[i+is][n];}
    fwrite(data,sizeof(float),ndata[0],p_binfile);
  }

#ifdef MHD
  /* Write out the interface magnetic fields */
  for (i=0;i<=ndata[0]; i++) {data[i] = (float)agrid->bx[i+is];}
  fwrite(data,sizeof(float),1+ndata[0],p_binfile);
#endif

  /* Close dump file, increment filename */

  fclose(p_binfile);

  /* Free the memory we malloc'd, but only once. */
  free(xgrid);
  data = NULL; 

#ifdef WRITE_DX_HEADER
  write_dx_header(agrid);
#endif /* WRITE_DX_HEADER */
  add1_2name(agrid->bin_file);
}

#endif /* ONE_D */

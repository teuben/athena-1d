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

/* Uncommenting the following define will cause the output routines to
   write the ghost cells as well as the computational cells.  This is
   usefull for debugging possible errors in setting the boundary
   conditions. */
/*
  #define WRITE_GHOST_CELLS
*/


#ifdef WRITE_DX_HEADER

#ifdef MHD
#define NUM_ARRAY (NVAR + 2)
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
  int nx, ny;
  char array_name[][18] = 
    {"Mass Density",
     "1-Momentum",
     "2-Momentum",
     "3-Momentum"
#ifdef ADIABATIC
     ,"Energy Density"
#endif
#ifdef MHD
     ,"1-Field","2-Field","3-Field"
     ,"1-Interface-Field","2-Interface-Field"
#endif
    };

  strcpy(FileName,agrid->bin_file);/* Copy the Data Dump File Name */
  strcat(FileName,".dx"); /* Append ".dx" for the header name */
#ifdef WRITE_GHOST_CELLS
  nx = agrid->ie1 - agrid->is1 + 9;
  ny = agrid->ie2 - agrid->is2 + 9;
#else
  nx = agrid->ie1 - agrid->is1 + 1;
  ny = agrid->ie2 - agrid->is2 + 1;
#endif
  nzones = nx*ny;

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

  fprintf(pfile,"# object 1 is the grid cell corner positions\n");
  fprintf(pfile,"object 1 class gridpositions counts %d %d\n",nx+1,ny+1);
  fprintf(pfile,"origin 0 0\n");
  fprintf(pfile,"delta  %g 0\n",agrid->dx);
  fprintf(pfile,"delta  0 %g\n\n",agrid->dy);

  fprintf(pfile,"# object 2 are the regular connections, quads\n");
  fprintf(pfile,"object 2 class gridconnections counts %d %d\n\n",nx+1,ny+1);

#ifdef MHD
  fprintf(pfile,"# object 3 is the cell x-interface positions\n");
#ifdef WRITE_GHOST_CELLS
  fprintf(pfile,"object 3 class gridpositions counts %d %d\n",nx,ny);
#else
  fprintf(pfile,"object 3 class gridpositions counts %d %d\n",nx+1,ny);
#endif
  fprintf(pfile,"origin 0 %g\n",0.5*agrid->dy);
  fprintf(pfile,"delta  %g 0\n",agrid->dx);
  fprintf(pfile,"delta  0 %g\n\n",agrid->dy);

  fprintf(pfile,"# object 4 is the cell y-interface positions\n");
  fprintf(pfile,"object 4 class gridpositions counts %d %d\n",nx,ny+1);
#ifdef WRITE_GHOST_CELLS
  fprintf(pfile,"object 4 class gridpositions counts %d %d\n",nx,ny);
#else
  fprintf(pfile,"object 4 class gridpositions counts %d %d\n",nx,ny+1);
#endif
  fprintf(pfile,"origin %g 0\n",0.5*agrid->dx);
  fprintf(pfile,"delta  %g 0\n",agrid->dx);
  fprintf(pfile,"delta  0 %g\n\n",agrid->dy);
#endif /* MHD */

  fprintf(pfile,"# objects 5 on are the data, which ");
  fprintf(pfile,"are in a one-to-one correspondence\n");
  fprintf(pfile,"# with the connections (\"dep\" on connections).\n\n");

  /* Initialize the byte offset to skip over the header info. */
  offset = 3*sizeof(int) + 2*sizeof(float);

  for(n=0; n<NVAR; n++){
    fprintf(pfile,"# %s\n",array_name[n]);
    fprintf(pfile,"object %d class array type float rank 0 items %d\n",
	    n+5,nzones);
    fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
    fprintf(pfile,"attribute \"dep\" string \"connections\"\n\n");
    offset += nzones*sizeof(float);
  }

#ifdef MHD
  fprintf(pfile,"# %s\n",array_name[NUM_ARRAY-2]);
#ifdef WRITE_GHOST_CELLS
  fprintf(pfile,"object %d class array type float rank 0 items %d\n",
	  n+5,nzones);
#else
  fprintf(pfile,"object %d class array type float rank 0 items %d\n",
	  n+5,nzones+ny);
#endif
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");

#ifndef WRITE_GHOST_CELLS
  /* If we are not writing out the ghost cells, there is one more
     column of interface fields. */
  offset += (nzones + ny)*sizeof(float);
#endif

  fprintf(pfile,"# %s\n",array_name[NUM_ARRAY-1]);
#ifdef WRITE_GHOST_CELLS
  fprintf(pfile,"object %d class array type float rank 0 items %d\n",
	  n+6,nzones);
#else
  fprintf(pfile,"object %d class array type float rank 0 items %d\n",
	  n+6,nzones+nx);
#endif
  fprintf(pfile,"data file \"%s\",%d\n",agrid->bin_file,offset);
  fprintf(pfile,"attribute \"dep\" string \"positions\"\n\n");
#endif /* MHD */

  for(n=0; n<72; n++) fputc((int)'#',pfile);

  fprintf(pfile,"\n\n# Next we construct a field for each of the objects above.\n");

  for(n=0; n<NVAR; n++){
    fprintf(pfile,"object \"%s\" class field\n",array_name[n]);
    fprintf(pfile,"component \"positions\" value 1\n");
    fprintf(pfile,"component \"connections\" value 2\n");
    fprintf(pfile,"component \"data\" value %d\n\n",n+5);
  }

#ifdef MHD
  fprintf(pfile,"object \"%s\" class field\n",array_name[NUM_ARRAY-2]);
  fprintf(pfile,"component \"positions\" value 3\n");
  fprintf(pfile,"component \"data\" value %d\n\n",NUM_ARRAY+3);

  fprintf(pfile,"object \"%s\" class field\n",array_name[NUM_ARRAY-1]);
  fprintf(pfile,"component \"positions\" value 4\n");
  fprintf(pfile,"component \"data\" value %d\n\n",NUM_ARRAY+4);
#endif /* MHD */

  for(n=0; n<72; n++) fputc((int)'#',pfile);

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



void binary_dump(struct grid_block *agrid){
  FILE *p_binfile;
  int ndata[3],is1,ie1,is2,ie2,i,j,n;
  float eos[2],*data;

  if((p_binfile = fopen(agrid->bin_file,"wb")) == NULL){
    fprintf(stderr,"[binary_dump]: Unable to open \"%s\" file\n",
	    agrid->bin_file);
    return;
  }
  is1 = agrid->is1;
  ie1 = agrid->ie1;
  is2 = agrid->is2;
  ie2 = agrid->ie2;
  /* Write number of zones and variables */

#ifdef WRITE_GHOST_CELLS
  ndata[0] = ie1-is1+9;
  ndata[1] = ie2-is2+9;
#else
  ndata[0] = ie1-is1+1;
  ndata[1] = ie2-is2+1;
#endif
  ndata[2] = NVAR;
  fwrite(ndata,sizeof(int),3,p_binfile);

  /* Write (gamma-1) and isothermal sound speed */

  eos[0] = (float)GAMM1 ;
  eos[1] = (float)ISOTHERMAL_C ;
  fwrite(eos,sizeof(float),2,p_binfile);

  /* Allocate Memory */

  if((data = (float *)malloc((1 + ndata[1])*sizeof(float))) == NULL){
    fprintf(stderr,"[binary_dump]: Unable to malloc memory to write binary file %s\n",agrid->bin_file);
    fclose(p_binfile);
    add1_2name(agrid->bin_file);
    return;
  }

  /* Write data */

  for(n=0; n<NVAR; n++){
    for(i=0; i<ndata[0]; i++){
      for(j=0; j<ndata[1]; j++){
#ifdef WRITE_GHOST_CELLS
	data[j] = (float)agrid->u[i+is1-4][j+is2-4][n];
#else
	data[j] = (float)agrid->u[i+is1][j+is2][n];
#endif
      }
      fwrite(data,sizeof(float),ndata[1],p_binfile);
    }
  }

#ifdef MHD
#ifdef WRITE_GHOST_CELLS
  /* Write out the interface magnetic fields */
  for(i=0; i<ndata[0]; i++){
    for(j=0; j<ndata[1]; j++){
      data[j] = (float)agrid->bx[i+is1-4][j+is2-4];
    }
    fwrite(data,sizeof(float),ndata[1],p_binfile);
  }

  for(i=0; i<ndata[0]; i++){
    for(j=0; j<ndata[1]; j++){
      data[j] = (float)agrid->by[i+is1-4][j+is2-4];
    }
    fwrite(data,sizeof(float),1+ndata[1],p_binfile);
  }
#else

  /* Write out the interface magnetic fields */
  for(i=0; i<=ndata[0]; i++){
    for(j=0; j<ndata[1]; j++){
      data[j] = (float)agrid->bx[i+is1][j+is2];
    }
    fwrite(data,sizeof(float),ndata[1],p_binfile);
  }

  for(i=0; i<ndata[0]; i++){
    for(j=0; j<=ndata[1]; j++){
      data[j] = (float)agrid->by[i+is1][j+is2];
    }
    fwrite(data,sizeof(float),1+ndata[1],p_binfile);
  }
#endif
#endif

  /* Close dump file, increment filename */
  fclose(p_binfile);

  /* Free the memory we malloc'd */
  free(data);

#ifdef WRITE_DX_HEADER
  write_dx_header(agrid);
#endif /* WRITE_DX_HEADER */
  add1_2name(agrid->bin_file);
}


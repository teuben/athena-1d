/*==============================================================================
 *                              Function PRINTD
 * Dumps scalar "history" variables in a formatted write
 *  scal[0] = time
 *  scal[1] = dt
 *  scal[2] = mass
 *  scal[3] = total energy
 *  scal[4] = 0.5*d*v1**2
 *  scal[5] = 0.5*d*v2**2
 *  scal[6] = 0.5*d*v3**2
 *  scal[7] = 0.5*b1**2
 *  scal[8] = 0.5*b2**2
 *  scal[9] = 0.5*b3**2
 * More variables can be added by increasing NSCAL=number of variables
 *============================================================================*/
#include <stdio.h>
#include "athena.def"
#include "athena.h"

void printd(struct grid_block *agrid)
{
#define NSCAL 10
  int is,ie,i=0;
  Real dvol, scal[NSCAL];
  FILE *p_hstfile;

   is = agrid->is;
   ie = agrid->ie;
   dvol = 1.0/(double)((ie-is+1));
   for (i=2; i<NSCAL; i++) {
      scal[i] = 0.0;
   }
 
/* Compute history variables */

   scal[0] = agrid->time;
   scal[1] = agrid->dt;
   for (i=is; i<=ie; i++) {
      scal[2] += agrid->u[i][0];
#ifdef ADIABATIC
      scal[3] += agrid->u[i][4];
#endif
      scal[4] += 0.5*SQR(agrid->u[i][1])/agrid->u[i][0];
      scal[5] += 0.5*SQR(agrid->u[i][2])/agrid->u[i][0];
      scal[6] += 0.5*SQR(agrid->u[i][3])/agrid->u[i][0];
#ifdef MHD
      scal[7] += 0.125*SQR(agrid->bx[i] + agrid->bx[i+1]);
      scal[8] += 0.25*SQR(agrid->u[i][NVAR-2]);
      scal[9] += 0.25*SQR(agrid->u[i][NVAR-1]);
#endif
   }
   for (i=2; i<NSCAL; i++) {
      scal[i] *= dvol;
   }

/* Write history file */

   p_hstfile = fopen(agrid->hst_file, "a");
   for (i=0; i<NSCAL; i++) {
      fprintf(p_hstfile,"% 14.6e",scal[i]);
   }
   fprintf(p_hstfile,"\n");
   fclose(p_hstfile);
}
#undef NSCAL

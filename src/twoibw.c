/*============================================================================*/
/*///////////////////////////// Function TWOIBW \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
/*                                                                            */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
void twoibw(FILE *p_input_file, struct grid_block *agrid)
{
/* Problem generator for two interacting blast waves (from WC)                */
/*                                                                            */
/*============================================================================*/
#include "prototypes.h"
char buf[120];
int i=0;
int is,ie;
/*============================================================================*/
/* Read idirect from 'athinput' */

/* setup dependent variables for shocktube in X-direction */

   is = agrid->is; ie = agrid->ie;

   for (i=is; i<=ie; i++) {
      agrid->u[i][0] = 1.0;
      agrid->u[i][1] = 0.0;
      agrid->u[i][2] = 0.0;
      agrid->u[i][3] = 0.0;
      if ((((double)(i-is) + 0.5)*agrid->dx) < 0.1) {
         agrid->u[i][4] = 1.0e3/GAMM1 ;
      } else if ((((double)(i-is) + 0.5)*agrid->dx) > 0.9) {
         agrid->u[i][4] = 1.0e2/GAMM1 ;
      } else {
         agrid->u[i][4] = 0.01/GAMM1 ;
      }
   }

}

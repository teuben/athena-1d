/*============================================================================*/
/*///////////////////////////// Function SHKSET \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
/*                                                                            */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

void shu_osher(FILE *p_input_file, struct grid_block *agrid)
{
/* Problem generator for 1-D Riemann problems                                 */
/*                                                                            */
/*============================================================================*/
  int i=0;
  int is,ie;
  Real dl,pl,ul,vl,wl,x[NX1];
/*============================================================================*/

   is = agrid->is; ie = agrid->ie;
   x[is-2] = -1.0 - 1.5*agrid->dx;
   for (i=is-1; i<=ie+2; i++) {
     x[i] = x[i-1]+agrid->dx;
   }

/* setup dependent variables for shocktube in X-direction */

   dl = 3.857143;
   ul = 2.629369;
   pl = 10.3333;
   wl = 0.0;
   vl = 0.0;
   for (i=is-2; i<=ie+2; i++) {
   if (x[i] < -0.8) {
      agrid->u[i][0] = dl;
      agrid->u[i][1] = ul*dl;
      agrid->u[i][2] = vl*dl;
      agrid->u[i][3] = wl*dl;
      agrid->u[i][4] = pl/GAMM1 
        + 0.5*dl*(ul*ul + vl*vl + wl*wl);
   } else {
      agrid->u[i][0] = 1.0 + 0.2*sin(5.0*PI*x[i]);
      agrid->u[i][1] = 0.0;
      agrid->u[i][2] = 0.0;
      agrid->u[i][3] = 0.0;
      agrid->u[i][4] = 1.0/GAMM1;
   }}
}

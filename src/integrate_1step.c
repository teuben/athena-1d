/*==============================================================================
 *                      Function INTEGRATE_1STEP 
 *
 * Updates the input "grid_block" structure pointed to by *agrid by one
 * timestep using Roe's scheme. The maximum eigenvalue (wave speed in the
 * linearized system) is found over all zones and
 * used to compute the next timestep.  This function also incre-
 * ments the "time" and "ncycles" stored in the grid_block
 *
 * Input: *agrid = pointer to "grid_block" structure
 *   Variables updated are: d,Sx,Sy,Sz,E,Bx,By,Bz,time,dt,ncycles
 *
 *============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
void integrate_1step(struct grid_block *agrid)
{
#include "prototypes.h"
int i=0,n,is,ie;
REAL dtodx;
REAL maxev_lr,maxev_roe,maxev_xswp=0.0,min_dtx=1.0e37;
REAL u[NXMAX][NVAR],wl[NXMAX][NVAR],wr[NXMAX][NVAR],f_u[NXMAX][NVAR];
REAL b_par[NXMAX+1];

   is = agrid->is;
   ie = agrid->ie;
   dtodx = agrid->dt/agrid->dx;

/* Extract 1-D pencil, packing variables into U[i] in correct order.  For */
/* X1-sweep: U = (d, Sx, Sy, Sz, E, By, Bz) */

   for (i=is-3; i<=ie+3; i++) {
      for (n=0; n<NVAR; n++) {
         u[i][n] = agrid->u[i][n];
      }
#ifdef MHD
      b_par[i] = agrid->bx[i];
#endif
   }
#ifdef MHD
   b_par[ie+4] = agrid->bx[ie+4];
#endif

/* Compute L and R states and fluxes, save maximum eigenvalue, and update */
/* hydro variables */

   maxev_lr  = lr_states (u,b_par,dtodx,is,ie,wl,wr    );
   maxev_roe = roe_fluxes(wl,wr,b_par,  is,ie,f_u);
   maxev_xswp = MAX(maxev_xswp,(MAX(maxev_lr,maxev_roe)));

   for (i=is; i<=ie; i++) {
   for (n=0; n<NVAR; n++) {
      agrid->u[i][n] -= dtodx*(f_u[i+1][n] - f_u[i][n]);
   }}

   min_dtx = agrid->dx/maxev_xswp;

/*================================= Finish up  ===============================*/

   agrid->time = agrid->time + agrid->dt;
   agrid->ncycles = agrid->ncycles + 1;
   agrid->dt = COURNO*min_dtx;
}

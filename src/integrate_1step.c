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
#include "prototypes.h"

#if defined ONE_D

void integrate_1step(struct grid_block *agrid)
{
  int i,n,is,ie;
  Real dtodx;
  Real maxev_lr,maxev_roe,maxev_xswp=0.0,min_dtx;
  Real u[NXMAX][NVAR],wl[NXMAX][NWAVE],wr[NXMAX][NWAVE],f_u[NXMAX][NWAVE];
  Real b_par[NXMAX];

  is = agrid->is;
  ie = agrid->ie;
  dtodx = agrid->dt/agrid->dx;

/* Extract 1-D pencil, packing variables into U[i] in correct order.  For */
/* X1-sweep: U = (d, Sx, Sy, Sz, E, Bx, By, Bz) */

  for (i=is-3; i<=ie+3; i++) {
    for (n=0; n<NVAR; n++) {
      u[i][n] = agrid->u[i][n];
    }
#ifdef MHD
    b_par[i] = agrid->bx[i];
#endif
  }

/* Compute L and R states and fluxes, save maximum eigenvalue, and update */
/* hydro variables */

  maxev_lr  = lr_states (u,b_par,dtodx,is,ie,wl,wr);
  maxev_roe = roe_fluxes(wl,wr,b_par,is,ie,f_u);
  maxev_xswp = MAX(maxev_xswp,(MAX(maxev_lr,maxev_roe)));

  for (i=is; i<=ie; i++) {
#ifdef MHD
    for (n=0; n<NVAR-3; n++) {
      agrid->u[i][n] -= dtodx*(f_u[i+1][n] - f_u[i][n]);
    }
    agrid->u[i][NVAR-2] -= dtodx*(f_u[i+1][NWAVE-2] - f_u[i][NWAVE-2]);
    agrid->u[i][NVAR-1] -= dtodx*(f_u[i+1][NWAVE-1] - f_u[i][NWAVE-1]);
#else
    for (n=0; n<NVAR; n++) {
      agrid->u[i][n] -= dtodx*(f_u[i+1][n] - f_u[i][n]);
    }
#endif
  }

  min_dtx = agrid->dx/maxev_xswp;

/*================================= Finish up  ===============================*/

  agrid->time = agrid->time + agrid->dt;
  agrid->ncycles = agrid->ncycles + 1;
  agrid->dt = COURNO*min_dtx;
}

#elif defined TWO_D

void integrate_1step(Grid_Block *agrid){
  int i,j,n,is1,ie1,is2,ie2;
  Real dtodx,dtody;
  Real maxev_lr,maxev_roe,maxev_xswp=0.0,min_dtx,maxev_yswp=0.0,min_dty;
  Real u[NXMAX][NVAR],wl[NXMAX][NWAVE],wr[NXMAX][NWAVE];
  Real b_par[NXMAX];
  Real f1_u[NXMAX][NXMAX][NWAVE]; /* NOTE: the order is f1_u[iy][ix][nvar] */
  Real f2_u[NXMAX][NXMAX][NWAVE]; /* NOTE: the order is f1_u[ix][iy][nvar] */
  Real tmp;

  is1 = agrid->is1;
  ie1 = agrid->ie1;
  dtodx = agrid->dt/agrid->dx;

  is2 = agrid->is2;
  ie2 = agrid->ie2;
  dtody = agrid->dt/agrid->dy;

/*================================= X1 Fluxes  ===============================*/

  /* Extract 1-D pencil, packing variables into U[i] in correct order. */
  /* For X1-sweep: U = (d, Sx, Sy, Sz, E, Bx, By, Bz) */
  /* f1_u(d, Sx, Sy, Sz, E, By, Bz) */
  for(j=is2; j<=ie2; j++){
    for(i=is1-3; i<=ie1+3; i++){
      for(n=0; n<NVAR; n++){
	u[i][n] = agrid->u[i][j][n];
      }
#ifdef MHD
      b_par[i] = agrid->bx[i][j];
#endif
    }

    /* Compute L and R states and fluxes and save maximum eigenvalue */

    maxev_lr  = lr_states(u,b_par,dtodx,is1,ie1,wl,wr);
    maxev_roe = roe_fluxes(wl,wr,b_par,is1,ie1,f1_u[j]);
    maxev_xswp = MAX(maxev_xswp,(MAX(maxev_lr,maxev_roe)));
  }
  min_dtx = agrid->dx/maxev_xswp;

/*================================= X2 Fluxes  ===============================*/

  /* Extract 1-D pencil, packing variables into U[i] in correct order. */
  /* For X2-sweep: U = (d, Sy, Sz, Sx, E, By, Bz, Bx) */
  /* f2_u(d, Sy, Sz, Sx, E, Bz, Bx) */
  for(i=is1; i<=ie1; i++){
    for(j=is2-3; j<=ie2+3; j++){
      u[j][0] = agrid->u[i][j][0];
      u[j][1] = agrid->u[i][j][2];
      u[j][2] = agrid->u[i][j][3];
      u[j][3] = agrid->u[i][j][1];
#if defined ADIABATIC
      u[j][4] = agrid->u[i][j][4];
#endif /* ADIABATIC */
#ifdef MHD
      u[j][NVAR-3] = agrid->u[i][j][NVAR-2];
      u[j][NVAR-2] = agrid->u[i][j][NVAR-1];
      u[j][NVAR-1] = agrid->u[i][j][NVAR-3];
      b_par[j] = agrid->by[i][j];
#endif
    }

    /* Compute L and R states and fluxes and save maximum eigenvalue */

    maxev_lr  = lr_states(u,b_par,dtody,is2,ie2,wl,wr);
    maxev_roe = roe_fluxes(wl,wr,b_par,is2,ie2,f2_u[i]);
    maxev_yswp = MAX(maxev_yswp,(MAX(maxev_lr,maxev_roe)));

    /* Permute the fluxes back to the natural order */
    /* After Perm. f2_u(d, Sx, Sy, Sz, E, Bx, Bz) */
    for(j=is2-3; j<=ie2+3; j++){
      /* Momenta */
      tmp = f2_u[i][j][3];
      f2_u[i][j][3] = f2_u[i][j][2];
      f2_u[i][j][2] = f2_u[i][j][1];
      f2_u[i][j][1] = tmp;
#ifdef MHD
      /* Electric fields */
      tmp = f2_u[i][j][NWAVE-1];
      f2_u[i][j][NWAVE-1] = f2_u[i][j][NWAVE-2];
      f2_u[i][j][NWAVE-2] = tmp;
#endif
    }
  }
  min_dty = agrid->dy/maxev_yswp;

/*=================================== Update =================================*/

  for(i=is1; i<=ie1; i++){
    for(j=is2; j<=ie2; j++){
#ifdef MHD
      for(n=0; n<NVAR-3; n++){
	agrid->u[i][j][n] -= dtodx*(f1_u[j][i+1][n] - f1_u[j][i][n])
	  + dtody*(f2_u[i][j+1][n] - f2_u[i][j][n]);
      }

      agrid->u[i][j][NVAR-3] -= 
	dtody*(f2_u[i][j+1][NWAVE-2] - f2_u[i][j][NWAVE-2]);

      agrid->u[i][j][NVAR-2] -= 
	dtodx*(f1_u[j][i+1][NWAVE-2] - f1_u[j][i][NWAVE-2]);

      agrid->u[i][j][NVAR-1] -= 
	dtodx*(f1_u[j][i+1][NWAVE-1] - f1_u[j][i][NWAVE-1])
	+ dtody*(f2_u[i][j+1][NWAVE-1] - f2_u[i][j][NWAVE-1]);

      /* Here's where we need some prescription for the interface
	 magnetic fields */
#else
      for(n=0; n<NVAR; n++){
	agrid->u[i][j][n] -= dtodx*(f1_u[j][i+1][n] - f1_u[j][i][n])
	  + dtody*(f2_u[i][j+1][n] - f2_u[i][j][n]);
      }
#endif
    }
  }

/*================================= Finish up  ===============================*/

  agrid->time = agrid->time + agrid->dt;
  agrid->ncycles = agrid->ncycles + 1;
  agrid->dt = COURNO*(min_dtx < min_dty ? min_dtx : min_dty);
}

#else /* THREE_D */

#endif

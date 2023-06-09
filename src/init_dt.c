/*==============================================================================
 *                            Function INIT_DT
 * Computes initial timestep for a new run using CFL condition.  Thereafter,
 * dt is computed from maximum eigenvalue computed during update.
 *============================================================================*/
#include <stdio.h>
#include <math.h>
#include "athena.def"
#include "athena.h"

#if defined ONE_D

Real init_dt(struct grid_block *agrid)
{
  int i;
  Real d,v1,v2,v3,vmax,bsq,va,qsq,p,cs,cmax=0.0,dt;
#ifdef MHD
  Real b1,b2,b3;
#endif /* MHD */

  for (i=agrid->is; i<=agrid->ie; i++) {
    d = agrid->u[i][0];
    v1 = agrid->u[i][1]/agrid->u[i][0];
    v2 = agrid->u[i][2]/agrid->u[i][0];
    v3 = agrid->u[i][3]/agrid->u[i][0];
    vmax = MAX(fabs(v1),MAX(fabs(v2),fabs(v3)));
#ifdef MHD
    b1 = agrid->bx[i];
    b2 = agrid->u[i][NVAR-2];
    b3 = agrid->u[i][NVAR-1];
    bsq = b1*b1 + b2*b2 + b3*b3;
    va = sqrt(bsq/d);
#else
    bsq = 0.0;
    va = 0.0;
#endif
#ifdef ADIABATIC
    qsq = v1*v1 + v2*v2 + v3*v3;
    p = MAX( GAMM1*(agrid->u[i][4] - 0.5*d*qsq - 0.5*bsq), TINY_NUMBER);
    cs = sqrt(GAMMA*p/d);
#else
    cs = ISOTHERMAL_C;
#endif
    cmax = MAX(cmax, cs + va + vmax);
  }

  printf("cmax = %e, dx = %e\n",cmax,agrid->dx);
  dt = COURNO*agrid->dx/cmax; 
  return(dt);
}

#elif defined TWO_D

Real init_dt(Grid_Block *agrid){
  int i,j;
  Real d,v1,v2,v3,vmax,bsq,va,qsq,p,cs,cmax=0.0,dt;
#ifdef MHD
  Real b1,b2,b3;
#endif /* MHD */

  for (i=agrid->is1; i<=agrid->ie1; i++) {
    for (j=agrid->is2; j<=agrid->ie2; j++) {
      d = agrid->u[i][j][0];
      v1 = agrid->u[i][j][1]/agrid->u[i][j][0];
      v2 = agrid->u[i][j][2]/agrid->u[i][j][0];
      v3 = agrid->u[i][j][3]/agrid->u[i][j][0];
      vmax = MAX(fabs(v1),MAX(fabs(v2),fabs(v3)));
#ifdef MHD
      b1 = agrid->bx[i][j];
      b2 = agrid->by[i][j];
      b3 = agrid->u[i][j][NVAR-1];
      bsq = b1*b1 + b2*b2 + b3*b3;
      va = sqrt(bsq/d);
#else
      bsq = 0.0;
      va = 0.0;
#endif
#ifdef ADIABATIC
      qsq = v1*v1 + v2*v2 + v3*v3;
      p = MAX( GAMM1*(agrid->u[i][j][4] - 0.5*d*qsq - 0.5*bsq), TINY_NUMBER);
      cs = sqrt(GAMMA*p/d);
#else
      cs = ISOTHERMAL_C;
#endif
      cmax = MAX(cmax, cs + va + vmax);
    }
  }

  printf("cmax = %e, dx = %e, dy = %e\n",cmax,agrid->dx,agrid->dy);
  dt = COURNO*(agrid->dx < agrid->dy ? agrid->dx : agrid->dy)/cmax; 

  return dt;
}

#else /* THREE_D */

#error : 3-D version of init_dt() not yet implemented

#endif

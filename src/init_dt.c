/*==============================================================================
 *                            Function INIT_DT
 * Computes initial timestep for a new run using CFL condition.  Thereafter,
 * dt is computed from maximum eigenvalue computed during update.
 *============================================================================*/
#include <stdio.h>
#include <math.h>
#include "athena.def"
#include "athena.h"

Real init_dt(struct grid_block *agrid)
{
  int i;
  Real d,v1,v2,v3,vmax,b1,b2,b3,bsq,va,qsq,p,cs,cmax=0.0,dt;

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

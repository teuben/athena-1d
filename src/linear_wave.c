/*============================================================================*/
/*/////////////////////////// Function FAST_WAVE \\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
/*                                                                            */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
REAL u0,amp;
int wave_flag;
void linear_wave(FILE *p_input_file, struct grid_block *agrid)
{
/* Problem generator for 1-D fast magnetosonic waves                          */
/*                                                                            */
/*============================================================================*/
#include "prototypes.h"
char buf[120];
int i=0;
int is,ie,n,m;
REAL d0,p0,v0,w0,bx0,by0,bz0,h0,xfact,yfact;
REAL x[NX1],ev[NVAR],rem[NVAR][NVAR],lem[NVAR][NVAR];
/*============================================================================*/

/* Read initial conditions */

   if (fscanf(p_input_file,"%i %[^\n]\n", &wave_flag, buf) != 2) {
      printf("Error reading wave_flag\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &amp, buf) != 2) {
      printf("Error reading amp\n"); exit(EXIT_FAILURE);
   }
   if (fscanf(p_input_file,"%lg %[^\n]\n", &u0, buf) != 2) {
      printf("Error reading u0\n"); exit(EXIT_FAILURE);
   }
   d0 = 1.0;
   p0 = 1.0/GAMMA;
   v0 = 0.0;
   w0 = 0.0;
   bx0 = 1.0;
   by0 = 1.0;
   bz0 = sqrt(1.25);
   for (n=0; n<NVAR; n++) {
   for (m=0; m<NVAR; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
   }}

/* Get eigenmatrix */

   xfact = 0.0;
   yfact = 1.0;
#if defined(ISOTHERMAL) AND defined(MHD)
   Eigensystem4_isothermal_mhd_InConsVars(d0,u0,v0,w0,
      bx0,by0,bz0,xfact,yfact,ev,rem,lem);
#endif
#if defined(ADIABATIC) AND defined(MHD)
   h0 = (p0/GAMM1+0.5*(bx0*bx0+by0*by0+bz0*bz0)+0.5*d0*(u0*u0+v0*v0+w0*w0))
	  + (p0+0.5*(bx0*bx0+by0*by0+bz0*bz0));
   Eigensystem4_adiabatic_mhd_InConsVars(d0,u0,v0,w0,h0,
      bx0,by0,bz0,xfact,yfact,ev,rem,lem);
#endif
   printf("Ux - Cf = %e\n",ev[0]);
   printf("Ux - Ca = %e\n",ev[1]);
   printf("Ux - Cs = %e\n",ev[2]);
   printf("Ux      = %e\n",ev[3]);
   printf("Ux + Cs = %e\n",ev[4]);
   printf("Ux + Ca = %e\n",ev[5]);
   printf("Ux + Cf = %e\n",ev[6]);

/* Initialize wave */

   is = agrid->is; ie = agrid->ie;
   x[is] = 0.5*(1.0/(double)(ie-is+1));
   for (i=is+1; i<=ie; i++) {x[i] = x[i-1] + 1.0/(double)(ie-is+1);}
   for (i=is; i<=ie; i++) {
      agrid->u[i][0] = d0 + amp*cos(2.0*PI*x[i])*rem[0][wave_flag];
      agrid->u[i][1] = d0*u0 + amp*cos(2.0*PI*x[i])*rem[1][wave_flag];
      agrid->u[i][2] = d0*v0 + amp*cos(2.0*PI*x[i])*rem[2][wave_flag];
      agrid->u[i][3] = d0*w0 + amp*cos(2.0*PI*x[i])*rem[3][wave_flag];
      agrid->bx[i] = bx0;
      agrid->u[i][NVAR-2] = by0 + amp*cos(2.0*PI*x[i])*rem[NVAR-2][wave_flag];
      agrid->u[i][NVAR-1] = bz0 + amp*cos(2.0*PI*x[i])*rem[NVAR-1][wave_flag];
#ifdef ADIABATIC
      agrid->u[i][4] = p0/GAMM1 
        + 0.5*(bx0*bx0 + by0*by0 + bz0*bz0)
        + 0.5*d0*(u0*u0 + v0*v0 + w0*w0)+amp*cos(2.0*PI*x[i])*rem[4][wave_flag];
#endif
   }
   agrid->bx[ie+1] = bx0;
}
REAL linear_wave_error(struct grid_block *agrid)
{
/* Computes L1-error in linear waves, ASSUMING WAVE HAS PROPAGATED AN INTEGER 
 * NUMBER OF PERIODS
 *
 *============================================================================*/
#include "prototypes.h"
int i=0;
int is,ie,n,m;
REAL d0,p0,v0,w0,bx0,by0,bz0,h0,xfact,yfact;
REAL x[NX1],ev[NVAR],rem[NVAR][NVAR],lem[NVAR][NVAR],error[NVAR];
/*============================================================================*/

printf("amp=%e, u0=%e, flag=%i\n",amp,u0,wave_flag);
   d0 = 1.0;
   p0 = 1.0/GAMMA;
   v0 = 0.0;
   w0 = 0.0;
   bx0 = 1.0;
   by0 = 1.0;
   bz0 = sqrt(1.25);
   for (n=0; n<NVAR; n++) {
   for (m=0; m<NVAR; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
   }}

/* Get eigenmatrix */

   xfact = 0.0;
   yfact = 1.0;
#if defined(ISOTHERMAL) AND defined(MHD)
   Eigensystem4_isothermal_mhd_InConsVars(d0,u0,v0,w0,
      bx0,by0,bz0,xfact,yfact,ev,rem,lem);
#endif
#if defined(ADIABATIC) AND defined(MHD)
   h0 = (p0/GAMM1+0.5*(bx0*bx0+by0*by0+bz0*bz0)+0.5*d0*(u0*u0+v0*v0+w0*w0))
	  + (p0+0.5*(bx0*bx0+by0*by0+bz0*bz0));
   Eigensystem4_adiabatic_mhd_InConsVars(d0,u0,v0,w0,h0,
      bx0,by0,bz0,xfact,yfact,ev,rem,lem);
#endif

   is = agrid->is; ie = agrid->ie;
   x[is] = 0.5*(1.0/(double)(ie-is+1));
   for (i=is+1; i<=ie; i++) {x[i] = x[i-1] + 1.0/(double)(ie-is+1);}
   for (n=0; n<NVAR; n++) {error[n] = 0.0;}
   for (i=is; i<=ie; i++) {
     error[0]+=fabs(agrid->u[i][0]-d0-amp*cos(2.0*PI*x[i])*rem[0][wave_flag]);
     error[1]+=fabs(agrid->u[i][1]-d0*u0-amp*cos(2.0*PI*x[i])*rem[1][wave_flag]);
     error[2]+=fabs(agrid->u[i][2]-d0*v0-amp*cos(2.0*PI*x[i])*rem[2][wave_flag]);
     error[3]+=fabs(agrid->u[i][3]-d0*w0-amp*cos(2.0*PI*x[i])*rem[3][wave_flag]);
     error[NVAR-2] += fabs(agrid->u[i][NVAR-2]
		     - by0 - amp*cos(2.0*PI*x[i])*rem[NVAR-2][wave_flag]);
     error[NVAR-1] += fabs(agrid->u[i][NVAR-1]
		     -bz0 - amp*cos(2.0*PI*x[i])*rem[NVAR-1][wave_flag]);
#ifdef ADIABATIC
     error[4] += fabs(agrid->u[i][4] - p0/GAMM1 
       - 0.5*(bx0*bx0 + by0*by0 + bz0*bz0)
       - 0.5*d0*(u0*u0 + v0*v0 + w0*w0)-amp*cos(2.0*PI*x[i])*rem[4][wave_flag]);
#endif
   }
   for (n=0; n<NVAR; n++) {error[n] = error[n]/(double)(ie-is+1);}
   printf("errors\n");
   printf("%i ",ie-is+1);
   for (n=0; n<NVAR; n++) {printf("%e ",error[n]);}
   printf("\n");
   for (n=0; n<NVAR; n++) {printf("rem[%i] = %e\n",n,rem[n][wave_flag]);}
}

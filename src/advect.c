/*============================================================================*/
/*////////////////////////////// Function ADVECT \\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/
/*                                                                            */
/* Problem generator for advection test.                                      */
/*                                                                            */
/*============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"

void advect(FILE *p_input_file, struct grid_block *agrid)
{
  char buf[120];
  int i=0,j=0,k=0,idirect;
  int is,ie,js,je,ks,ke;
  Real x[NX1],y[NX2],z[NX3];
  Real radx,rady,radz,rad2,amp,da,ua,va,wa,bxa,bya,bza;
/*============================================================================*/
/* Read idirect from 'athinput' */

   if (fscanf(p_input_file,"%i %[^\n]\n", &idirect, buf) != 2) {
      printf("Error reading idirect\n"); exit(EXIT_FAILURE);
   }

/* setup grid centers, starting at (0,0,0) */

   is = agrid->is; ie = agrid->ie;
   js = agrid->js; je = agrid->je;
   ks = agrid->ks; ke = agrid->ke;
   for(i=is;i<=ie;i++) x[i]=((float)(i-is)+0.5)*agrid->dx;
   for(j=js;j<=je;j++) y[j]=((float)(j-js)+0.5)*agrid->dy;
   for(k=ks;k<=ke;k++) z[k]=((float)(k-ks)+0.5)*agrid->dz;

/* Setup Gaussian centered at R=0.3, and square pulse for X=[0.7,0.9] */

   amp=1.0e-1;
   ua = 0.0;
   va = 0.0;
   wa = 0.0;
   bxa = 0.0;
   bya = 0.0;
   bza = 0.0;
   if (idirect == 1) {ua=1.0;}
   if (idirect == 2) {va=1.0;}
   if (idirect == 3) {wa=1.0;}
   if (idirect == 12) {ua=-1.0; va=-1.0;}
   for (k=ks; k<=ke; k++) {
   for (j=js; j<=je; j++) {
   for (i=is; i<=ie; i++) {
      radx = x[i]-0.3;
      rady = y[j]-0.3;
      radz = z[k]-0.3;
      if (idirect == 1) {
        agrid->d [k][j][i] = 1.0 + amp*exp(radx*radx/(-0.01));
        if (radx < 0.6 AND radx > 0.4) {agrid->d [k][j][i] = 1.0 + amp;}
      }
      if (idirect == 2) {
        agrid->d [k][j][i] = 1.0 + amp*exp(rady*rady/(-0.01));
        if (rady < 0.6 AND rady > 0.4) {agrid->d [k][j][i] = 1.0 + amp;}
      }
      if (idirect == 3) {
        agrid->d [k][j][i] = 1.0 + amp*exp(radz*radz/(-0.01));
        if (radz < 0.6 AND radz > 0.4) {agrid->d [k][j][i] = 1.0 + amp;}
      }
      if (idirect == 12) {
        rad2 = radx*radx + rady*rady;
        agrid->d [k][j][i] = 1.0 + amp*exp(rad2/(-0.01));
        if (radx < 0.6 AND radx > 0.4 AND rady < 0.6 AND rady > 0.4)
           {agrid->d [k][j][i] = 1.0 + amp;}
      }
      agrid->sx[k][j][i] = agrid->d [k][j][i]*ua;
      agrid->sy[k][j][i] = agrid->d [k][j][i]*va;
      agrid->sz[k][j][i] = agrid->d [k][j][i]*wa;
#ifdef MHD
      agrid->bx[k][j][i] = bxa;
      agrid->by[k][j][i] = bya;
      agrid->bz[k][j][i] = bza;
      if (i == ie) agrid->bx[k][j][i+1] = bxa;
      if (j == je) agrid->by[k][j+1][i] = bya;
      if (k == ke) agrid->bz[k+1][j][i] = bza;
#endif
#ifdef ADIABATIC
      agrid->e [k][j][i] = 1.0e-5/GAMM1 
#ifdef MHD
        + 0.5*(bxa*bxa + bya*bya + bza*bza)
#endif
        + 0.5*agrid->d [k][j][i]*(ua*ua + va*va + wa*wa);
#endif
   }}}
}

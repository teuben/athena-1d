/*==============================================================================
 *                            Function ROE_FLUXES
 * Computes fluxes of conserved quantities along a 1-D slice using Roe's
 * linearization (1981, JCP, 43, 357).  When Roe's linearization fails because
 * of negative density or total energy in the intermediate states, then the
 * Harten-Lax-vanLeer-Einfeldt (HLLE) solver of Einfeldt et al. (1991, JCP, 92,
 * 273) is used.
 *
 * Input Arguments:
 *    wl, wr = L/R-states of PRIMITIVE variables in 1-D slice at cell faces
 *    b1 = component of B in direction of slice at cell faces
 *    ibegin,iend = starting and ending indices of zone centers in slice
 * wl, wr, and b1 must be defined over [ibegin:iend+1]
 *
 * Output Arguments:
 *    f = fluxes of CONSERVED variables at cell edges over [ibegin:iend+1]
 *
 * Returns:
 *   maxevroe = maximum absolute value of eigenvalues (effective wave speeds) of
 *   Roe matrix over every zone face along slice (used to compute timestep)
 *
 *============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

REAL roe_fluxes(REAL wl[NXMAX][NVAR], REAL wr[NXMAX][NVAR],
   REAL b1[NXMAX+1], int ibegin, int iend, REAL f[NXMAX][NVAR])
{
  REAL sqrtdl,sqrtdr,sdlpdr;
  REAL droe,v1roe,v2roe,v3roe,b2roe,b3roe,hroe,x,y,pbl=0.0,pbr=0.0,el,er;
  REAL ev[NVAR],rem[NVAR][NVAR],lem[NVAR][NVAR],maxevroe=0.0;
  REAL ul[NVAR],ur[NVAR],du[NVAR],a[NVAR],u_inter[NVAR],p_inter=0.0;
  REAL fl[NVAR],fr[NVAR];
  int i,n,m,hlle_flag;
  REAL asq,vaxsq=0.0,qsq,cfsq,cfl,cfr,bp,bm;

   for (n=0; n<NVAR; n++) {
   for (m=0; m<NVAR; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
   }}

/* Begin loop over all zone faces */
   for (i=ibegin; i<=iend+1; i++) {

/*--- Step 1. ------------------------------------------------------------------
 * Compute Roe-averaged data from left- and right-states
 */
      sqrtdl = sqrt(wl[i][0]);
      sqrtdr = sqrt(wr[i][0]);
      sdlpdr = sqrtdl + sqrtdr;

      droe  = sqrtdl*sqrtdr;
      v1roe = (sqrtdl*wl[i][1] + sqrtdr*wr[i][1])/sdlpdr;
      v2roe = (sqrtdl*wl[i][2] + sqrtdr*wr[i][2])/sdlpdr;
      v3roe = (sqrtdl*wl[i][3] + sqrtdr*wr[i][3])/sdlpdr;
/*
 * The Roe average of the magnetic field is defined with sqrtdl[r] reversed
 * compared to the Roe average of the other variables.  The numerical
 * factors X and Y are needed to compute the eigenvectors (eqs. B15,B16)
 */
#ifdef MHD
      pbl = 0.5*(SQR(wl[i][NVAR-2]) + SQR(wl[i][NVAR-1]) + b1[i]*b1[i]);
      pbr = 0.5*(SQR(wr[i][NVAR-2]) + SQR(wr[i][NVAR-1]) + b1[i]*b1[i]);
      b2roe = (sqrtdr*wl[i][NVAR-2] + sqrtdl*wr[i][NVAR-2])/sdlpdr;
      b3roe = (sqrtdr*wl[i][NVAR-1] + sqrtdl*wr[i][NVAR-1])/sdlpdr;
      x = 0.5*((b2roe*b2roe - wl[i][NVAR-2]*wr[i][NVAR-2])
              +(b3roe*b3roe - wl[i][NVAR-1]*wr[i][NVAR-1]))/droe;
      y = 0.5*(wl[i][0] + wr[i][0])/droe;
#endif
#ifdef ADIABATIC
/*
 * Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
 * rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
 */
      el = wl[i][4]/GAMM1 + pbl;
      el += 0.5*wl[i][0]*(SQR(wl[i][1]) + SQR(wl[i][2]) + SQR(wl[i][3]));
      er = wr[i][4]/GAMM1 + pbr;
      er += 0.5*wr[i][0]*(SQR(wr[i][1]) + SQR(wr[i][2]) + SQR(wr[i][3]));
      hroe  = ((el+wl[i][4]+pbl)/sqrtdl + (er+wr[i][4]+pbr)/sqrtdr)/sdlpdr;
#endif

/*--- Step 2. ------------------------------------------------------------------
 * Compute eigenvalues and eigenmatrices using Roe-averaged values
 */
#if defined(ISOTHERMAL) AND defined(HYDRO)
      Eigensystem4_isothermal_hydro_InConsVars(v1roe,v2roe,v3roe,ev,rem,lem);
#endif
#if defined(ADIABATIC) AND defined(HYDRO)
      Eigensystem4_adiabatic_hydro_InConsVars(v1roe,v2roe,v3roe,hroe,ev,rem,lem);
#endif
#if defined(ISOTHERMAL) AND defined(MHD)
      Eigensystem4_isothermal_mhd_InConsVars(droe,v1roe,v2roe,v3roe,
	 b1[i],b2roe,b3roe,x,y,ev,rem,lem);
#endif
#if defined(ADIABATIC) AND defined(MHD)
      Eigensystem4_adiabatic_mhd_InConsVars(droe,v1roe,v2roe,v3roe,hroe,
	 b1[i],b2roe,b3roe,x,y,ev,rem,lem);
#endif
      maxevroe = MAX(maxevroe,(MAX(fabs(ev[0]),fabs(ev[NVAR-1]))));

/*--- Step 3. ------------------------------------------------------------------
 * Check that the density and pressure in the intermediate states are positive.
 * If not, set hlle_flag=1 if d_inter<0; hlle_flag=2 if p_inter<0
 */
      hlle_flag = 0;
/* Compute L- and R-states of conserved variables */
      ul[0] = wl[i][0];
      ul[1] = wl[i][0]*wl[i][1];
      ul[2] = wl[i][0]*wl[i][2];
      ul[3] = wl[i][0]*wl[i][3];
      ur[0] = wr[i][0];
      ur[1] = wr[i][0]*wr[i][1];
      ur[2] = wr[i][0]*wr[i][2];
      ur[3] = wr[i][0]*wr[i][3];
#ifdef ADIABATIC
      ul[4] = el;
      ur[4] = er;
#endif
#ifdef MHD
      ul[NVAR-2] = wl[i][NVAR-2];
      ul[NVAR-1] = wl[i][NVAR-1];
      ur[NVAR-2] = wr[i][NVAR-2];
      ur[NVAR-1] = wr[i][NVAR-1];
#endif
/* Now compute intermediate states of conserved variables using Roe solver */
      for (n=0; n<NVAR; n++) {du[n] = ur[n] - ul[n];}
      for (n=0; n<NVAR; n++) {
	 a[n] = 0.0;
         for (m=0; m<NVAR; m++) {a[n] += lem[n][m]*du[m];}
      }
      for (n=0; n<NVAR; n++) {u_inter[n] = ul[n];}
      for (n=0; n<NVAR-1; n++) {
         for (m=0; m<NVAR; m++) {u_inter[m] += a[n]*rem[m][n];}
         if (u_inter[0] < 0.0) {hlle_flag=1;}
#ifdef ADIABATIC
	 p_inter = u_inter[4] - 0.5*
	       (SQR(u_inter[1])+SQR(u_inter[2])+SQR(u_inter[3]))/u_inter[0];
#ifdef MHD
	 p_inter -= 0.5*(SQR(u_inter[NVAR-2])+SQR(u_inter[NVAR-1])+SQR(b1[i]));
#endif
	 if (p_inter < 0.0) {hlle_flag=2;}
#endif
      }

/*--- Step 4. ------------------------------------------------------------------
 * Compute fluxes using L/R states 
 */
      fl[0] = ul[1];
      fr[0] = ur[1];
      fl[1] = ul[1]*wl[i][1];
      fr[1] = ur[1]*wr[i][1];
      fl[2] = ul[1]*wl[i][2];
      fr[2] = ur[1]*wr[i][2];
      fl[3] = ul[1]*wl[i][3];
      fr[3] = ur[1]*wr[i][3];

#ifdef ADIABATIC
      fl[1] += wl[i][4];
      fr[1] += wr[i][4];
      fl[4] = (ul[4] + wl[i][4])*wl[i][1];
      fr[4] = (ur[4] + wr[i][4])*wr[i][1];
#else
      fl[1] += wl[i][0]*(ISOTHERMAL_C_SQ);
      fr[1] += wr[i][0]*(ISOTHERMAL_C_SQ);
#endif

#ifdef MHD
      fl[1] -= 0.5*(b1[i]*b1[i] - SQR(wl[i][NVAR-2]) - SQR(wl[i][NVAR-1]));
      fr[1] -= 0.5*(b1[i]*b1[i] - SQR(wr[i][NVAR-2]) - SQR(wr[i][NVAR-1]));
      fl[2] -= b1[i]*wl[i][NVAR-2];
      fr[2] -= b1[i]*wr[i][NVAR-2];
      fl[3] -= b1[i]*wl[i][NVAR-1];
      fr[3] -= b1[i]*wr[i][NVAR-1];
#ifdef ADIABATIC
      fl[4] += (pbl*wl[i][1] - b1[i]*(b1[i]*wl[i][1] + wl[i][NVAR-2]*wl[i][2]
				                    + wl[i][NVAR-1]*wl[i][3]));
      fr[4] += (pbr*wr[i][1] - b1[i]*(b1[i]*wr[i][1] + wr[i][NVAR-2]*wr[i][2]
                                                    + wr[i][NVAR-1]*wr[i][3]));
#endif
      fl[NVAR-2] = wl[i][NVAR-2]*wl[i][1] - b1[i]*wl[i][2];
      fr[NVAR-2] = wr[i][NVAR-2]*wr[i][1] - b1[i]*wr[i][2];
      fl[NVAR-1] = wl[i][NVAR-1]*wl[i][1] - b1[i]*wl[i][3];
      fr[NVAR-1] = wr[i][NVAR-1]*wr[i][1] - b1[i]*wr[i][3];
#endif

/*--- Step 5. ------------------------------------------------------------------
 * Compute fluxes at interface using either the Roe or the HLLE solver
 */
      if (hlle_flag != 0) {

/* Use HLLE solver */

         printf("hlle_flag=%i, HLLE solver used for i=%i\n",hlle_flag,i); 
#ifdef ADIABATIC
         asq = GAMMA*wl[i][4]/wl[i][0]; 
#else
         asq = ISOTHERMAL_C_SQ ;
#endif
#ifdef MHD
         vaxsq = b1[i]*b1[i]/wl[i][0];
#endif
         qsq = asq + pbl/wl[i][0];
         cfsq = 0.5*(qsq + sqrt(qsq*qsq-4.*asq*vaxsq));
         cfl = sqrt(cfsq);
#ifdef ADIABATIC
         asq = GAMMA*wr[i][4]/wr[i][0]; 
#else
         asq = ISOTHERMAL_C_SQ ;
#endif
#ifdef MHD
         vaxsq = b1[i]*b1[i]/wr[i][0];
#endif
         qsq = asq + pbr/wr[i][0];
         cfsq = 0.5*(qsq + sqrt(qsq*qsq-4.*asq*vaxsq));
         cfr = sqrt(cfsq);
         bp = MAX(MAX(ev[NVAR-1],(wr[i][1] + cfr)), 0.0);
         bm = MIN(MIN(ev[0]     ,(wl[i][1] - cfl)), 0.0);
         for (n=0; n<NVAR; n++) {
            f[i][n] = ((bp*fl[n]-bm*fr[n]) + bp*bm*(ur[n]-ul[n]))/(bp-bm);
         }

/* Use Roe solver */

      } else {
         for (n=0; n<NVAR; n++) {
         f[i][n] = 0.5*(fl[n] + fr[n]);
         for (m=0; m<NVAR; m++) {
            f[i][n] -= 0.5*fabs(ev[m])*a[m]*rem[n][m];
         }}
      }

/* end loop over i */

   }
   return(maxevroe);
}

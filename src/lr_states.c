/*==============================================================================
 *                           Function LR_STATES
 * Computes left- and right-states of PRIMITIVE variables using the Piecewise
 * Parabolic Method of Colella & Woodward (1985, JCP,54,174) along a 1-D slice
 *
 * Input Arguments:
 *   u = CONSERVED variables at cell centers along 1-D slice
 *   b1 = component of b in direction of slice at cell faces
 *   dtodx = dt/dx
 *   ib,ie = beginning and ending indices of zone centers in slice
 * Including ghost zones, U must be initialized over [ib-3:ie+3], and b1 over
 * [ib-3:ie+4]
 *
 * Output Arguments:
 *   wl,wr= L/R-states of PRIMITIVE variables at cell faces over [ib:ie+1]
 *
 * Returns:
 *   maxevlr = maximum absolute value of eigenvalues over [ib:ie] (for dt)
 *
 *============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

Real lr_states(Real u[NXMAX][NVAR], Real b1[NXMAX], Real dtodx,
	       int ib, int ie, Real wl[NXMAX][NWAVE], Real wr[NXMAX][NWAVE])
{
  int i,n,m;
  Real w[NXMAX][NWAVE],pb=0.0,ke;
  Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE],b1c=0.0,maxevlr=0.0;
  Real ev_ip1[NWAVE],rem_ip1[NWAVE][NWAVE],lem_ip1[NWAVE][NWAVE];
  Real dwc[NWAVE],dwl[NWAVE],dwr[NWAVE],dwg[NWAVE],dwm[NXMAX][NWAVE];
  Real dac[NWAVE],dal[NWAVE],dar[NWAVE],dag[NWAVE],da[NWAVE];
#ifdef THIRD_ORDER
  Real wim1h[NXMAX][NWAVE];
#endif
  Real wlv[NWAVE],wrv[NWAVE],dw[NWAVE],w6[NWAVE];
  Real qa,qd,qe,qx;

  for (n=0; n<NWAVE; n++) {
    for (m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
      rem_ip1[n][m] = 0.0;
      lem_ip1[n][m] = 0.0;
    }
  }

/*--- Step 1. ------------------------------------------------------------------
 * Transform to primitive variables, W=(d,v1,v2,v3,[P],[b2,b3])
 */
   for (i=ib-3; i<=ie+3; i++) {
      w[i][0] = u[i][0];
      w[i][1] = u[i][1]/u[i][0];
      w[i][2] = u[i][2]/u[i][0];
      w[i][3] = u[i][3]/u[i][0];
#ifdef MHD
      w[i][NWAVE-2] = u[i][NVAR-2];
      w[i][NWAVE-1] = u[i][NVAR-1];
      pb = 0.5*(u[i][NVAR-2]*u[i][NVAR-2] + u[i][NVAR-1]*u[i][NVAR-1] + 
		u[i][NVAR-3]*u[i][NVAR-3] + 
		ONE_3RD*(b1[i] - u[i][NVAR-3])*(b1[i] - u[i][NVAR-3]) );
#endif
#ifdef ADIABATIC
      ke = 0.5*w[i][0]*(w[i][1]*w[i][1] + w[i][2]*w[i][2] + w[i][3]*w[i][3]);
      w[i][4] = MAX(GAMM1*(u[i][4] - ke - pb),TINY_NUMBER);
#endif
   }

/*--- Step 2. ------------------------------------------------------------------
 * Compute TVD linear slopes at ib-2 and ib-1.  The right eigenvectors are
 * stored as COLUMNS of rem; the left eigenvectors are stored as ROWS of lem
 */

   for (i=ib-2; i<=ib-1; i++) {
#if defined(ISOTHERMAL) AND defined(HYDRO)
      Eigensystem4_isothermal_hydro_InPrimVars(w[i][0],w[i][1],w[i][2],
	 w[i][3],ev,rem,lem);
#endif
#if defined(ADIABATIC) AND defined(HYDRO)
      Eigensystem4_adiabatic_hydro_InPrimVars(w[i][0],w[i][1],w[i][2],
	 w[i][3],w[i][4],ev,rem,lem);
#endif
#if defined(ISOTHERMAL) AND defined(MHD)
      b1c = u[i][NVAR-3];
      Eigensystem4_isothermal_mhd_InPrimVars(w[i][0],w[i][1],w[i][2],
         w[i][3],b1c,w[i][4],w[i][5],ev,rem,lem);
#endif
#if defined(ADIABATIC) AND defined(MHD)
      b1c = u[i][NVAR-3];
      Eigensystem4_adiabatic_mhd_InPrimVars(w[i][0],w[i][1],w[i][2],
         w[i][3],w[i][4],b1c,w[i][5],w[i][6],ev,rem,lem);
#endif
      maxevlr = MAX(maxevlr,(MAX(fabs(ev[0]),fabs(ev[NWAVE-1]))));

/* Compute centered, L/R, and van Leer differences of primitive variables */

      for (n=0; n<NWAVE; n++) {
         dwc[n] = w[i+1][n] - w[i-1][n];
         dwl[n] = w[i  ][n] - w[i-1][n];
         dwr[n] = w[i+1][n] - w[i  ][n];
         if (dwl[n]*dwr[n] > 0.0) {
            dwg[n] = 2.0*dwl[n]*dwr[n]/(dwl[n]+dwr[n]);
         } else {
            dwg[n] = 0.0;
         }
      }

/* Project differences of primitive variables along characteristics */

      for (n=0; n<NWAVE; n++) {
         dac[n] = lem[n][0]*dwc[0];
         dal[n] = lem[n][0]*dwl[0];
         dar[n] = lem[n][0]*dwr[0];
         dag[n] = lem[n][0]*dwg[0];
         for (m=1; m<NWAVE; m++) {
            dac[n] += lem[n][m]*dwc[m];
            dal[n] += lem[n][m]*dwl[m];
            dar[n] += lem[n][m]*dwr[m];
            dag[n] += lem[n][m]*dwg[m];
         }
      }

/* Apply monotonicity constraint to projections along characteristics */

      for (n=0; n<NWAVE; n++) {
         da[n] = 0.0;
         if (dal[n]*dar[n] > 0.0) {
           da[n] = SIGN(dac[n])*MIN(2.0*MIN(    fabs(dal[n]),fabs(dar[n])),
                                        MIN(0.5*fabs(dac[n]),fabs(dag[n])) );
         }
      }

/* Compute monotonic slopes of primitive variables from projections */

      for (n=0; n<NWAVE; n++) {
         dwm[i][n] = da[0]*rem[n][0];
         for (m=1; m<NWAVE; m++) {
            dwm[i][n] += da[m]*rem[n][m];
         }
      }
   }

/*--- Step 3. ------------------------------------------------------------------
 * For third order scheme, construct interpolant in primitive variables at edge
 * of cell ib-1 ("W[ib-1/2]", CW eqn 1.6) using linear TVD slopes at ib-2 and
 * ib-1 computed in Step 2.
 */

#ifdef THIRD_ORDER
   for (n=0; n<NWAVE; n++) {
     wim1h[ib-1][n] = .5*(w[ib-1][n]+w[ib-2][n])-(dwm[ib-1][n]-dwm[ib-2][n])/6.;
   }
#endif

/*--- Step 4. ------------------------------------------------------------------
 * Compute L/R-states of primitive variables over entire 1D slice. This
 * requires repeating steps 2 and 3 above to get wim1h for ib:ie+1, and then
 * integrating the interpolation polynomial to get L/R-states at cell faces.
 *
 * At the start of the loop, rem and lem still store values at i=ib-1 computed
 * at the end of Step 2.  For each i, the eigensystem at i+1 is stored in
 * rem_ip1 and lem_ip1.  At the end of the loop rem[lem] is then set to
 * rem_ip1[lem_ip1] in preparation for the next iteration.
 */

   for (i=ib-1; i<=ie+1; i++) {
#if defined(ISOTHERMAL) AND defined(HYDRO)
      Eigensystem4_isothermal_hydro_InPrimVars(w[i+1][0],w[i+1][1],
	 w[i+1][2],w[i+1][3],ev_ip1,rem_ip1,lem_ip1);
#endif
#if defined(ADIABATIC) AND defined(HYDRO)
      Eigensystem4_adiabatic_hydro_InPrimVars(w[i+1][0],w[i+1][1],
	 w[i+1][2],w[i+1][3],w[i+1][4],ev_ip1,rem_ip1,lem_ip1);
#endif
#if defined(ISOTHERMAL) AND defined(MHD)
      b1c = u[i+1][NVAR-3];
      Eigensystem4_isothermal_mhd_InPrimVars(w[i+1][0],w[i+1][1],
	 w[i+1][2],w[i+1][3],b1c,w[i+1][4],w[i+1][5],ev_ip1,rem_ip1,lem_ip1);
#endif
#if defined(ADIABATIC) AND defined(MHD)
      b1c = u[i+1][NVAR-3];
      Eigensystem4_adiabatic_mhd_InPrimVars(w[i+1][0],w[i+1][1],w[i+1][2]
         ,w[i+1][3],w[i+1][4],b1c,w[i+1][5],w[i+1][6],ev_ip1,rem_ip1,lem_ip1);
#endif
      maxevlr = MAX(maxevlr,(MAX(fabs(ev_ip1[0]),fabs(ev_ip1[NWAVE-1]))));

/* Compute centered, L/R, and van Leer differences of primitive variables */

      for (n=0; n<NWAVE; n++) {
         dwc[n] = w[i+2][n] - w[i  ][n];
         dwl[n] = w[i+1][n] - w[i  ][n];
         dwr[n] = w[i+2][n] - w[i+1][n];
         if (dwl[n]*dwr[n] > 0.0) {
            dwg[n] = 2.0*dwl[n]*dwr[n]/(dwl[n]+dwr[n]);
         } else {
            dwg[n] = 0.0;
         }
      }

/* Project differences of primitive variables along characteristics */

      for (n=0; n<NWAVE; n++) {
         dac[n] = lem_ip1[n][0]*dwc[0];
         dal[n] = lem_ip1[n][0]*dwl[0];
         dar[n] = lem_ip1[n][0]*dwr[0];
         dag[n] = lem_ip1[n][0]*dwg[0];
         for (m=1; m<NWAVE; m++) {
            dac[n] += lem_ip1[n][m]*dwc[m];
            dal[n] += lem_ip1[n][m]*dwl[m];
            dar[n] += lem_ip1[n][m]*dwr[m];
            dag[n] += lem_ip1[n][m]*dwg[m];
         }
      }

/* Apply monotonicity constraint to projections along characteristics */

      for (n=0; n<NWAVE; n++) {
         da[n] = 0.0;
         if (dal[n]*dar[n] > 0.0) {
            da[n] = SIGN(dac[n])*MIN(2.0*MIN(    fabs(dal[n]),fabs(dar[n])),
                                         MIN(0.5*fabs(dac[n]),fabs(dag[n])) );
         }
      }

/* Compute monotonic slopes of primitive variables from projections */

      for (n=0; n<NWAVE; n++) {
         dwm[i+1][n] = da[0]*rem_ip1[n][0];
         for (m=1; m<NWAVE; m++) {
            dwm[i+1][n] += da[m]*rem_ip1[n][m];
         }
      }
/*
 * For third order scheme, construct interpolant in primitive variables at cell
 * edges (wim1h="W[i minus one-half]", CW eqn 1.6).  Then set values at left- 
 * and right-edges of cell i (the a_L,i and a_R,i).
 */
#ifdef THIRD_ORDER
      for (n=0; n<NWAVE; n++) {
         wim1h[i+1][n] = 0.5*(w[i+1][n]+w[i][n]) - (dwm[i+1][n]-dwm[i][n])/6.0;
      }
      for (n=0; n<NWAVE; n++) {
         wlv[n] = wim1h[i  ][n];
         wrv[n] = wim1h[i+1][n];
      }
#endif
#ifdef SECOND_ORDER
      for (n=0; n<NWAVE; n++) {
         wlv[n] = w[i][n] - 0.5*dwm[i][n];
         wrv[n] = w[i][n] + 0.5*dwm[i][n];
      }
#endif
/* Ensure the a_L,i and a_R,i lie between cell-center values */

      for (n=0; n<NWAVE; n++) {
         wlv[n] = MAX(MIN(w[i][n],w[i-1][n]),wlv[n]);
         wlv[n] = MIN(MAX(w[i][n],w[i-1][n]),wlv[n]);
         wrv[n] = MAX(MIN(w[i][n],w[i+1][n]),wrv[n]);
         wrv[n] = MIN(MAX(w[i][n],w[i+1][n]),wrv[n]);
      }

/* Monotonize again (CW eqn 1.10) */

      for (n=0; n<NWAVE; n++) {
         qa = (wrv[n]-w[i][n])*(w[i][n]-wlv[n]);
         qd = wrv[n]-wlv[n];
         qe = 6.0*(w[i][n] - 0.5*(wlv[n]+wrv[n]));
         if (qa <= 0.0) {
            wlv[n] = w[i][n];
            wrv[n] = w[i][n];
         } else if (qd*(qd - qe) < 0.0) {
            wlv[n] = 3.0*w[i][n] - 2.0*wrv[n];
         } else if (qd*(qd + qe) < 0.0) {
            wrv[n] = 3.0*w[i][n] - 2.0*wlv[n];
         }

/* Compute coefficients of interpolation parabolae (CW eqn 1.5) */

         dw[n] = wrv[n] - wlv[n];
         w6[n] = 6.0*(w[i][n] - 0.5*(wlv[n] + wrv[n]));
      }
/*
 * Compute left(right)-state by integrating interpolation parabolae over
 * domain of dependence defined by max(min) eigenvalue (CW eqn 1.12).  Note
 * left(right)-state at interface i+1/2 (i-1/2) is on the right(left)-side of
 * cell i
 */
      qx = TWO_3RDS*MAX(ev[NWAVE-1],0.0)*dtodx;
      for (n=0; n<NWAVE; n++) {
         wl[i+1][n] = wrv[n] - 0.75*qx*(dw[n] - (1.0 - qx)*w6[n]);
      }
      qx = -TWO_3RDS*MIN(ev[0],0.0)*dtodx;
      for (n=0; n<NWAVE; n++) {
         wr[i][n] = wlv[n] + 0.75*qx*(dw[n] + (1.0 - qx)*w6[n]);
      }
/*
 * Then subtract amount of each wave m that does not reach the interface
 * during timestep (CW eqn 3.5ff)
 */
      for (n=0; n<NWAVE-1; n++) {
      if (ev[n] > 0.) {
         qa  = 0.0;
         for (m=0; m<NWAVE; m++) {
            qa += lem[n][m]*0.5*dtodx*
             ( (ev[NWAVE-1]-ev[n])*(dw[m]-w6[m]) + dtodx*TWO_3RDS*w6[m]*
               (ev[NWAVE-1]*ev[NWAVE-1] - ev[n]*ev[n]) );
         }
         for (m=0; m<NWAVE; m++) {wl[i+1][m] += qa*rem[m][n];}
      }}

      for (n=1; n<NWAVE; n++) {
      if (ev[n] < 0.) {
         qa = 0.0;
         for (m=0; m<NWAVE; m++) {
            qa += lem[n][m]*0.5*dtodx*
             ( (ev[0]-ev[n])*(dw[m]+w6[m]) + dtodx*TWO_3RDS*w6[m]*
               (ev[0]*ev[0] - ev[n]*ev[n]) );
         }
         for (m=0; m<NWAVE; m++) {wr[i][m] += qa*rem[m][n];}
      }}

/* Save eigenvalues and eigenmatrices at i+1 for use in next iteration */

      for (m=0; m<NWAVE; m++) {
      ev[m] = ev_ip1[m];
      for (n=0; n<NWAVE; n++) {
         rem[m][n] = rem_ip1[m][n];
         lem[m][n] = lem_ip1[m][n];
      }}
   }

/*--- End of Step 4 (Loop over i). ---------------------------------------------
 * Return maximum eigenvalue and exit
 */

   return(maxevlr);
}

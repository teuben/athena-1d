/*============================================================================*
 *                            Function LR_FLUXES
 * Computes fluxes of CONSERVED variables from left- and right-states
 *
 * Input Arguments:
 *   ul,ur = L/R-states of CONSERVED variables
 *   wl,wr = L/R-states of PRIMITIVE variables
 * Both ul/ur and wl/wr are used to simplify expressions and save flops
 *   b1 = longitudinal component of B
 *   pbl,pbr = magnetic pressure in L/R states
 *
 * Output Arguments:
 *   fl,fr = fluxes of CONSERVED variables computed from L/R states
 *
 *============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
#include "prototypes.h"

void lr_fluxes(REAL ul[NVAR], REAL ur[NVAR], REAL wl[NVAR], REAL wr[NVAR],
   REAL b1, REAL pbl, REAL pbr, REAL fl[NVAR], REAL fr[NVAR])
{
   
   fl[0] = ul[1];
   fr[0] = ur[1];
   fl[1] = ul[1]*wl[1];
   fr[1] = ur[1]*wr[1];
   fl[2] = ul[1]*wl[2];
   fr[2] = ur[1]*wr[2];
   fl[3] = ul[1]*wl[3];
   fr[3] = ur[1]*wr[3];

#ifdef HYDRO
#ifdef ADIABATIC
   fl[1] += wl[4];
   fr[1] += wr[4];
   fl[4] = (ul[4] + wl[4])*wl[1];
   fr[4] = (ur[4] + wr[4])*wr[1];
#else
   fl[1] += wl[0]*(ISOTHERMAL_C_SQ);
   fr[1] += wr[0]*(ISOTHERMAL_C_SQ);
#endif
#endif

#ifdef MHD
   fl[1] -= 0.5*(b1*b1 - wl[NVAR-2]*wl[NVAR-2] - wl[NVAR-3]*wl[NVAR-3]);
   fr[1] -= 0.5*(b1*b1 - wr[NVAR-2]*wr[NVAR-2] - wr[NVAR-3]*wr[NVAR-3]);
   fl[2] -= b1*wl[NVAR-2];
   fr[2] -= b1*wr[NVAR-2];
   fl[3] -= b1*wl[NVAR-3];
   fr[3] -= b1*wr[NVAR-3];
#ifdef ADIABATIC
   fl[4] += pbl*wl[1] - b1*(b1*wl[1] + wl[NVAR-2]*wl[2] + wl[NVAR-3]*wl[3]);
   fr[4] += pbr*wr[1] - b1*(b1*wr[1] + wr[NVAR-2]*wr[2] + wr[NVAR-3]*wr[3]);
#endif
   fl[NVAR-2] = b2l[NVAR-2]*wl[1] - b1*wl[2];
   fr[NVAR-2] = b2r[NVAR-2]*wr[1] - b1*wr[2];
   fl[NVAR-1] = b3l[NVAR-3]*wl[1] - b1*wl[3];
   fr[NVAR-1] = b3r[NVAR-3]*wr[1] - b1*wr[3];
#endif
}

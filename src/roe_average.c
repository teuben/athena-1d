/*============================================================================*/
 *                            Function ROE_AVERAGE
 * Computes Roe averages of left- and right-states
 *
 * Input Arguments:
 *   wl,wr = L/R-states of PRIMITIVE variables
 *   b1 = longitudinal component of B
 *
 * Output Arguments:
 *   droe,v1roe,...,b3roe = Roe-averaged states
 *   x,y = numerical factors needed to compute eigenvalues
 *   pbl,pbr = magnetic pressure in L/R states (needed to compute fluxes)
 *   el,er = total energy in L/R states (needed to compute fluxes)
 *
 *============================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena.def"
#include "athena.h"
void roe_average(REAL wl[NVAR], REAL wr[NVAR], REAL b1, REAL droe,
   REAL v1roe, REAL v2roe, REAL v3roe, REAL hroe, REAL b2roe, REAL b3roe
   REAL x, REAL y, REAL pbl, REAL pbr, REAL el, REAL er)
{
#include "prototypes.h"
REAL sqrtdl,sqrtdr,sdlpdr;
   sqrtdl = sqrt(wl[0]);
   sqrtdr = sqrt(wr[0]);
   sdlpdr = sqrtdl + sqrtdr;

   droe  = sqrtdl*sqrtdr;
   v1roe = (sqrtdl*wl[1] + sqrtdr*wr[1])/sdlpdr;
   v2roe = (sqrtdl*wl[2] + sqrtdr*wr[2])/sdlpdr;
   v3roe = (sqrtdl*wl[3] + sqrtdr*wr[3])/sdlpdr;
/*
 * The Roe average of the magnetic field is defined with sqrtdl[r] reversed
 * compared to the Roe average of the other variables.  The numerical
 * factors X and Y are needed to compute the eigenvectors (eqs. B15,B16)
 */
#ifdef MHD
   pbl = 0.5*(wl[NVAR-2]*wl[NVAR-2] + wl[NVAR-1]*wl[NVAR-1] + b1*b1)
   pbr = 0.5*(wr[NVAR-2]*wr[NVAR-2] + wr[NVAR-1]*wr[NVAR-1] + b1*b1)
   b2roe = (sqrtdr*wl[NVAR-2] + sqrtdl*wr[NVAR-2])/sdlpdr;
   b3roe = (sqrtdr*wl[NVAR-1] + sqrtdl*wr[NVAR-1])/sdlpdr;
   x = 0.5*((b2roe*b2roe - wl[NVAR-2]*wr[NVAR-2])
           +(b3roe*b3roe - wl[NVAR-1]*wr[NVAR-1]))/droe;
   y = 0.5*(wl[0] + wr[0])/droe;
#endif
#ifdef ADIABATIC
/*
 * Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
 * rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
 */
   el = wl[4]/GAMM1 + 0.5*wl[0]*(wl[1]*wl[1]+wl[2]*wl[2]+wl[3]*wl[3]) + pbl;
   er = wr[4]/GAMM1 + 0.5*wr[0]*(wr[1]*wr[1]+wr[2]*wr[2]+wr[3]*wr[3]) + pbr;
   hroe  = ((el+wl[4])/sqrtdl + (er+wr[4])/sqrtdr)/sdlpdr;
#endif
}

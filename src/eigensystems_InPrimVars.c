/*==============================================================================
 *           File eigensystems_InPrimVars (Multiple Functions)
 * This file contains functions to evaluate the eigenvalues, and left- and
 * right-eigenvectors of "Roe's matrix A" for the linearized system in the 
 * PRIMITIVE variables, i.e. W,t = AW,x, where W=(d,vx,vy,vz,[P],[By,Bz]).
 * The eigenvalues are returned through the argument list as a vector of length
 * NVAR.  The eigenvectors are returned as matrices of size (NVAR)x(NVAR), with
 * right-eigenvectors stored as COLUMNS (so R_i = right_eigenmatrix[*][i]), and
 * left-eigenvectors stored as ROWS (so L_i = left_eigenmatrix[i][*]).
 *
 * To improve performance components of the eigenvectors which are zero
 * are not set here (eigenmatrices must be initialized to zero in calling
 * routine).   However, for completeness statements which set these values
 * are included, but are commented out.
 *
 *============================================================================*/
#include <math.h>
#include <float.h>
#include "athena.def"
#include "athena.h"
/*----------------------------- ISOTHERMAL HYDRO -------------------------------
 *
 * Input: v1,v2,v3 = components of velocity
 * Output: eigenvalues[NVAR]
 *         right_eigenmatrix[NVAR,NVAR], left_eigenmatrix[NVAR,NVAR];
 */
void Eigensystem4_isothermal_hydro_InPrimVars(Real d, Real v1, Real v2, Real v3,
   Real eigenvalues[NVAR],
   Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR])
{
/* Compute eigenvalues (eq. A5) */

   eigenvalues[0] = v1 - ISOTHERMAL_C;
   eigenvalues[1] = v1;
   eigenvalues[2] = v1;
   eigenvalues[3] = v1 + ISOTHERMAL_C;

/* Right-eigenvectors, stored as COLUMNS (eq. A3 sans 4th row and 4th column) */

   right_eigenmatrix[0][0] = 1.0;
   right_eigenmatrix[1][0] = -(ISOTHERMAL_C)/d;
/* right_eigenmatrix[2][0] = 0.0; */
   right_eigenmatrix[3][0] = ISOTHERMAL_C_SQ;

   right_eigenmatrix[0][1] = 1.0;
/* right_eigenmatrix[1][1] = 0.0; */
/* right_eigenmatrix[2][1] = 0.0; */
/* right_eigenmatrix[3][1] = 0.0; */

/* right_eigenmatrix[0][2] = 0.0; */
/* right_eigenmatrix[1][2] = 0.0; */
   right_eigenmatrix[2][2] = 1.0;
/* right_eigenmatrix[3][2] = 0.0; */

   right_eigenmatrix[0][3] = 1.0;
   right_eigenmatrix[1][3] = (ISOTHERMAL_C)/d;
/* right_eigenmatrix[2][3] = 0.0; */
   right_eigenmatrix[3][3] = ISOTHERMAL_C_SQ;

/* Left-eigenvectors, stored as ROWS (eq. A4 sans 4th row and 4th column) */

/* left_eigenmatrix[0][0] = 0.0; */
   left_eigenmatrix[0][1] = -0.5*d/(ISOTHERMAL_C);
/* left_eigenmatrix[0][2] = 0.0; */
   left_eigenmatrix[0][3] = 0.5/(ISOTHERMAL_C_SQ);

   left_eigenmatrix[1][0] = 1.0;
/* left_eigenmatrix[1][1] = 0.0; */
/* left_eigenmatrix[1][2] = 0.0; */
   left_eigenmatrix[1][3] = -1.0/(ISOTHERMAL_C_SQ);

/* left_eigenmatrix[2][0] = 0.0; */
/* left_eigenmatrix[2][1] = 0.0; */
   left_eigenmatrix[2][2] = 1.0;
/* left_eigenmatrix[2][3] = 0.0; */

/* left_eigenmatrix[3][0] = 0.0; */
   left_eigenmatrix[3][1] = 0.5*d/(ISOTHERMAL_C);
/* left_eigenmatrix[3][2] = 0.0; */
   left_eigenmatrix[3][3] = 0.5/(ISOTHERMAL_C_SQ);

}
/*---------------------------- ADIABATIC HYDRO ---------------------------------
 *
 * Input: d,v1,v2,v3,p = primitive variables
 * Output: eigenvalues[NVAR]
 *         right_eigenmatrix[NVAR,NVAR], left_eigenmatrix[NVAR,NVAR];
 */
void Eigensystem4_adiabatic_hydro_InPrimVars(Real d, Real v1, Real v2, Real v3,
   Real p, Real eigenvalues[NVAR],
   Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR])
{
  Real asq,a;
  asq = GAMMA*p/d;
  a = sqrt(asq);

/* Compute eigenvalues (eq. A2) */

   eigenvalues[0] = v1 - a;
   eigenvalues[1] = v1;
   eigenvalues[2] = v1;
   eigenvalues[3] = v1;
   eigenvalues[4] = v1 + a;

/* Right-eigenvectors, stored as COLUMNS (eq. A3) */

   right_eigenmatrix[0][0] = 1.0;
   right_eigenmatrix[1][0] = -a/d;
/* right_eigenmatrix[2][0] = 0.0; */
/* right_eigenmatrix[3][0] = 0.0; */
   right_eigenmatrix[4][0] = asq;

   right_eigenmatrix[0][1] = 1.0;
/* right_eigenmatrix[1][1] = 0.0; */
/* right_eigenmatrix[2][1] = 0.0; */
/* right_eigenmatrix[3][1] = 0.0; */
/* right_eigenmatrix[4][1] = 0.0; */

/* right_eigenmatrix[0][2] = 0.0; */
/* right_eigenmatrix[1][2] = 0.0; */
   right_eigenmatrix[2][2] = 1.0;
/* right_eigenmatrix[3][2] = 0.0; */
/* right_eigenmatrix[4][2] = 0.0; */

/* right_eigenmatrix[0][3] = 0.0; */
/* right_eigenmatrix[1][3] = 0.0; */
/* right_eigenmatrix[2][3] = 0.0; */
   right_eigenmatrix[3][3] = 1.0;
/* right_eigenmatrix[4][3] = 0.0; */

   right_eigenmatrix[0][4] = 1.0;
   right_eigenmatrix[1][4] = -right_eigenmatrix[1][0];
/* right_eigenmatrix[2][4] = 0.0; */
/* right_eigenmatrix[3][4] = 0.0; */
   right_eigenmatrix[4][4] = asq;

/* Left-eigenvectors, stored as ROWS (eq. A4) */

/* left_eigenmatrix[0][0] = 0.0; */
   left_eigenmatrix[0][1] = -0.5*d/a;
/* left_eigenmatrix[0][2] = 0.0; */
/* left_eigenmatrix[0][3] = 0.0; */
   left_eigenmatrix[0][4] = 0.5/asq;

   left_eigenmatrix[1][0] = 1.0;
/* left_eigenmatrix[1][1] = 0.0; */
/* left_eigenmatrix[1][2] = 0.0; */
/* left_eigenmatrix[1][3] = 0.0; */
   left_eigenmatrix[1][4] = -1.0/asq;

/* left_eigenmatrix[2][0] = 0.0; */
/* left_eigenmatrix[2][1] = 0.0; */
   left_eigenmatrix[2][2] = 1.0;
/* left_eigenmatrix[2][3] = 0.0; */
/* left_eigenmatrix[2][4] = 0.0; */

/* left_eigenmatrix[3][0] = 0.0; */
/* left_eigenmatrix[3][1] = 0.0; */
/* left_eigenmatrix[3][2] = 0.0; */
   left_eigenmatrix[3][3] = 1.0;
/* left_eigenmatrix[3][4] = 0.0; */

/* left_eigenmatrix[4][0] = 0.0; */
   left_eigenmatrix[4][1] = -left_eigenmatrix[0][1];
/* left_eigenmatrix[4][2] = 0.0; */
/* left_eigenmatrix[4][3] = 0.0; */
   left_eigenmatrix[4][4] = left_eigenmatrix[0][4];
}
/*----------------------------- ISOTHERMAL MHD ---------------------------------
 *
 * Input: d,v1,v2,v3,b1,b2,b3 = density, velocities, and B field
 * Output: eigenvalues[NVAR]
 *         right_eigenmatrix[NVAR,NVAR], left_eigenmatrix[NVAR,NVAR];
 */
void Eigensystem4_isothermal_mhd_InPrimVars(Real d, Real v1, Real v2, Real v3,
   Real b1, Real b2, Real b3, Real eigenvalues[NVAR],
   Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR])
{
  Real btsq,vaxsq,q_starsq,cfsq,cf,cssq,cs,bt,bet2,bet3,alpha_f,alpha_s;
  Real sqrtd,s,qf,qs,af,as,vax,nf,ns,af_prime,as_prime;
  btsq  = b2*b2 + b3*b3;

/* Compute fast- and slow-magnetosonic speeds (eq. A8 and A9) */

   vaxsq = b1*b1/d;
   q_starsq  = (ISOTHERMAL_C_SQ) + vaxsq + btsq/d;
   cfsq = 0.5*(q_starsq + sqrt(q_starsq*q_starsq-4.*(ISOTHERMAL_C_SQ)*vaxsq));
   cf  = sqrt(cfsq);
   cssq = 0.5*(q_starsq - sqrt(q_starsq*q_starsq-4.*(ISOTHERMAL_C_SQ)*vaxsq));
   if (cssq <= 0.0) {
      cssq = 0.0;
   }
   cs  = sqrt(cssq);

/* Compute beta(s) (eq A13) */

   bt  = sqrt(btsq);
   if (bt == 0.0) {
      bet2 = ONE_OVER_SQRT2 ;
      bet3 = ONE_OVER_SQRT2 ;
   } else {
      bet2 = b2/bt;
      bet3 = b3/bt;
   }

/* Compute alpha(s) (eq A12) */

   if ((cfsq-cssq) == 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
   } else if ( (ISOTHERMAL_C_SQ - cssq) <= 0.0) {
      alpha_f = 0.0;
      alpha_s = 1.0;
   } else if ( (cfsq - ISOTHERMAL_C_SQ) <= 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
   } else {
      alpha_f = sqrt((ISOTHERMAL_C_SQ - cssq)/(cfsq - cssq));
      alpha_s = sqrt((cfsq - ISOTHERMAL_C_SQ)/(cfsq - cssq));
   }

/* Compute Q(s) and A(s) (eq. A11), etc. */

   sqrtd = sqrt(d);
   s = SIGN(b1);
   qf = cf*alpha_f*s;
   qs = cs*alpha_s*s;
   af = (ISOTHERMAL_C)*alpha_f*sqrtd;
   as = (ISOTHERMAL_C)*alpha_s*sqrtd;

/* Compute eigenvalues (eq. A17) */

   vax = sqrt(vaxsq);
   eigenvalues[0] = v1 - cf;
   eigenvalues[1] = v1 - vax;
   eigenvalues[2] = v1 - cs;
   eigenvalues[3] = v1 + cs;
   eigenvalues[4] = v1 + vax;
   eigenvalues[5] = v1 + cf;

/* Right-eigenvectors, stored as COLUMNS (eq A10 sans 5th row and 4th column) */

   right_eigenmatrix[0][0] = d*alpha_f;
   right_eigenmatrix[1][0] = -cf*alpha_f;
   right_eigenmatrix[2][0] = qs*bet2;
   right_eigenmatrix[3][0] = qs*bet3;
   right_eigenmatrix[4][0] = as*bet2;
   right_eigenmatrix[5][0] = as*bet3;

/* right_eigenmatrix[0][1] = 0.0; */
/* right_eigenmatrix[1][1] = 0.0; */
   right_eigenmatrix[2][1] = -bet3;
   right_eigenmatrix[3][1] = bet2;
   right_eigenmatrix[4][1] = -bet3*s*sqrtd;
   right_eigenmatrix[5][1] = bet2*s*sqrtd;

   right_eigenmatrix[0][2] = d*alpha_s;
   right_eigenmatrix[1][2] = -cs*alpha_s;
   right_eigenmatrix[2][2] = -qf*bet2;
   right_eigenmatrix[3][2] = -qf*bet3;
   right_eigenmatrix[4][2] = -af*bet2;
   right_eigenmatrix[5][2] = -af*bet3;

   right_eigenmatrix[0][3] = right_eigenmatrix[0][2];
   right_eigenmatrix[1][3] = -right_eigenmatrix[1][2];
   right_eigenmatrix[2][3] = -right_eigenmatrix[2][2];
   right_eigenmatrix[3][3] = -right_eigenmatrix[3][2];
   right_eigenmatrix[4][3] = right_eigenmatrix[4][2];
   right_eigenmatrix[5][3] = right_eigenmatrix[5][2];

/* right_eigenmatrix[0][4] = 0.0; */
/* right_eigenmatrix[1][4] = 0.0; */
   right_eigenmatrix[2][4] = bet3;
   right_eigenmatrix[3][4] = -bet2;
   right_eigenmatrix[4][4] = right_eigenmatrix[4][1];
   right_eigenmatrix[5][4] = right_eigenmatrix[5][1];

   right_eigenmatrix[0][5] = right_eigenmatrix[0][0];
   right_eigenmatrix[1][5] = -right_eigenmatrix[1][0];
   right_eigenmatrix[2][5] = -right_eigenmatrix[2][0];
   right_eigenmatrix[3][5] = -right_eigenmatrix[3][0];
   right_eigenmatrix[4][5] = right_eigenmatrix[4][0];
   right_eigenmatrix[5][5] = right_eigenmatrix[5][0];

/* Left-eigenvectors, stored as ROWS (eq A14 sans 4th row and 5th column) */

   nf = 1.0/((ISOTHERMAL_C_SQ)*(1.0+alpha_s*alpha_s));
   ns = 1.0/((ISOTHERMAL_C_SQ)*(1.0+alpha_f*alpha_f));
   qf = ns*qf;
   qs = nf*qs;
   af_prime = ns*af/d;
   as_prime = nf*as/d;

/* left_eigenmatrix[0][0] = 0.0; */
   left_eigenmatrix[0][1] = -nf*cf*alpha_f;
   left_eigenmatrix[0][2] = qs*bet2;
   left_eigenmatrix[0][3] = qs*bet3;
   left_eigenmatrix[0][4] = as_prime*bet2;
   left_eigenmatrix[0][5] = as_prime*bet3;

/* left_eigenmatrix[1][0] = 0.0; */
/* left_eigenmatrix[1][1] = 0.0; */
   left_eigenmatrix[1][2] = -0.5*bet3;
   left_eigenmatrix[1][3] = 0.5*bet2;
   left_eigenmatrix[1][4] = -0.5*bet3*s/sqrtd;
   left_eigenmatrix[1][5] = 0.5*bet2*s/sqrtd;

/* left_eigenmatrix[2][0] = 0.0; */
   left_eigenmatrix[2][1] = -ns*cs*alpha_s;
   left_eigenmatrix[2][2] = -qf*bet2;
   left_eigenmatrix[2][3] = -qf*bet3;
   left_eigenmatrix[2][4] = -af_prime*bet2;
   left_eigenmatrix[2][5] = -af_prime*bet3;

/* left_eigenmatrix[3][0] = 0.0; */
   left_eigenmatrix[3][1] = -left_eigenmatrix[2][1];
   left_eigenmatrix[3][2] = -left_eigenmatrix[2][2];
   left_eigenmatrix[3][3] = -left_eigenmatrix[2][3];
   left_eigenmatrix[3][4] = left_eigenmatrix[2][4];
   left_eigenmatrix[3][5] = left_eigenmatrix[2][5];

/* left_eigenmatrix[4][0] = 0.0; */
/* left_eigenmatrix[4][1] = 0.0; */
   left_eigenmatrix[4][2] = -left_eigenmatrix[1][2];
   left_eigenmatrix[4][3] = -left_eigenmatrix[1][3];
   left_eigenmatrix[4][4] = left_eigenmatrix[1][4];
   left_eigenmatrix[4][5] = left_eigenmatrix[1][5];

/* left_eigenmatrix[5][0] = 0.0; */
   left_eigenmatrix[5][1] = -left_eigenmatrix[0][1];
   left_eigenmatrix[5][2] = -left_eigenmatrix[0][2];
   left_eigenmatrix[5][3] = -left_eigenmatrix[0][3];
   left_eigenmatrix[5][4] = left_eigenmatrix[0][4];
   left_eigenmatrix[5][5] = left_eigenmatrix[0][5];
}
/*-------------------------------- ADIABATIC MHD -------------------------------
 *
 * Input: d,v1,v2,v3,p,b1,b2,b3 = density, velocities, pressure, and B field
 * Output: eigenvalues[NVAR]
 *         right_eigenmatrix[NVAR,NVAR], left_eigenmatrix[NVAR,NVAR];
 */
void Eigensystem4_adiabatic_mhd_InPrimVars(Real d, Real v1, Real v2, Real v3,
   Real p, Real b1, Real b2, Real b3, Real eigenvalues[NVAR],
   Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR])
{
  Real btsq,asq,vaxsq,q_starsq,cfsq,cf,cssq,cs,bt,bet2,bet3,alpha_f,alpha_s;
  Real sqrtd,s,a,qf,qs,af,as,vax,na,af_prime,as_prime;
  btsq  = b2*b2 + b3*b3;

/* Compute fast- and slow-magnetosonic speeds (eq. A8 and A9) */

   asq = GAMMA*p/d;
   vaxsq = b1*b1/d;
   q_starsq  = asq + (vaxsq + btsq/d);
   cfsq = 0.5*(q_starsq + sqrt(q_starsq*q_starsq - 4.0*asq*vaxsq));
   cf  = sqrt(cfsq);
   cssq = 0.5*(q_starsq - sqrt(q_starsq*q_starsq - 4.0*asq*vaxsq));
   if (cssq <= 0.0) {
      cssq = 0.0;
   }
   cs  = sqrt(cssq);

/* Compute beta(s) (eq A13) */

   bt  = sqrt(btsq);
   if (bt == 0.0) {
      bet2 = ONE_OVER_SQRT2 ;
      bet3 = ONE_OVER_SQRT2 ;
   } else {
      bet2 = b2/bt;
      bet3 = b3/bt;
   }

/* Compute alpha(s) (eq A12) */

   if ((cfsq-cssq) == 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
   } else if ( (asq - cssq) <= 0.0) {
      alpha_f = 0.0;
      alpha_s = 1.0;
   } else if ( (cfsq - asq) <= 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
   } else {
      alpha_f = sqrt((asq - cssq)/(cfsq - cssq));
      alpha_s = sqrt((cfsq - asq)/(cfsq - cssq));
   }

/* Compute Q(s) and A(s) (eq. A11), etc. */

   sqrtd = sqrt(d);
   s = SIGN(b1);
   a = sqrt(asq);
   qf = cf*alpha_f*s;
   qs = cs*alpha_s*s;
   af = a*alpha_f*sqrtd;
   as = a*alpha_s*sqrtd;

/* Compute eigenvalues (eq. A7) */

   vax = sqrt(vaxsq);
   eigenvalues[0] = v1 - cf;
   eigenvalues[1] = v1 - vax;
   eigenvalues[2] = v1 - cs;
   eigenvalues[3] = v1;
   eigenvalues[4] = v1 + cs;
   eigenvalues[5] = v1 + vax;
   eigenvalues[6] = v1 + cf;

/* Right-eigenvectors, stored as COLUMNS (eq. A10) */

   right_eigenmatrix[0][0] = d*alpha_f;
   right_eigenmatrix[1][0] = -cf*alpha_f;
   right_eigenmatrix[2][0] = qs*bet2;
   right_eigenmatrix[3][0] = qs*bet3;
   right_eigenmatrix[4][0] = d*asq*alpha_f;
   right_eigenmatrix[5][0] = as*bet2;
   right_eigenmatrix[6][0] = as*bet3;

/* right_eigenmatrix[0][1] = 0.0; */
/* right_eigenmatrix[1][1] = 0.0; */
   right_eigenmatrix[2][1] = -bet3;
   right_eigenmatrix[3][1] = bet2;
/* right_eigenmatrix[4][1] = 0.0; */
   right_eigenmatrix[5][1] = -bet3*s*sqrtd;
   right_eigenmatrix[6][1] = bet2*s*sqrtd;

   right_eigenmatrix[0][2] = d*alpha_s;
   right_eigenmatrix[1][2] = -cs*alpha_s;
   right_eigenmatrix[2][2] = -qf*bet2;
   right_eigenmatrix[3][2] = -qf*bet3;
   right_eigenmatrix[4][2] = d*asq*alpha_s;
   right_eigenmatrix[5][2] = -af*bet2;
   right_eigenmatrix[6][2] = -af*bet3;

   right_eigenmatrix[0][3] = 1.0;
/* right_eigenmatrix[1][3] = 0.0; */
/* right_eigenmatrix[2][3] = 0.0; */
/* right_eigenmatrix[3][3] = 0.0; */
/* right_eigenmatrix[4][3] = 0.0; */
/* right_eigenmatrix[5][3] = 0.0; */
/* right_eigenmatrix[6][3] = 0.0; */

   right_eigenmatrix[0][4] = right_eigenmatrix[0][2];
   right_eigenmatrix[1][4] = -right_eigenmatrix[1][2];
   right_eigenmatrix[2][4] = -right_eigenmatrix[2][2];
   right_eigenmatrix[3][4] = -right_eigenmatrix[3][2];
   right_eigenmatrix[4][4] = right_eigenmatrix[4][2];
   right_eigenmatrix[5][4] = right_eigenmatrix[5][2];
   right_eigenmatrix[6][4] = right_eigenmatrix[6][2];

/* right_eigenmatrix[0][5] = 0.0; */
/* right_eigenmatrix[1][5] = 0.0; */
   right_eigenmatrix[2][5] = bet3;
   right_eigenmatrix[3][5] = -bet2;
/* right_eigenmatrix[4][5] = 0.0; */
   right_eigenmatrix[5][5] = right_eigenmatrix[5][1];
   right_eigenmatrix[6][5] = right_eigenmatrix[6][1];

   right_eigenmatrix[0][6] = right_eigenmatrix[0][0];
   right_eigenmatrix[1][6] = -right_eigenmatrix[1][0];
   right_eigenmatrix[2][6] = -right_eigenmatrix[2][0];
   right_eigenmatrix[3][6] = -right_eigenmatrix[3][0];
   right_eigenmatrix[4][6] = right_eigenmatrix[4][0];
   right_eigenmatrix[5][6] = right_eigenmatrix[5][0];
   right_eigenmatrix[6][6] = right_eigenmatrix[6][0];

/* Left-eigenvectors, stored as ROWS (eq. A14) */

   na = 0.5/asq;
   qf = na*qf;
   qs = na*qs;
   af_prime = na*af/d;
   as_prime = na*as/d;

/* left_eigenmatrix[0][0] = 0.0; */
   left_eigenmatrix[0][1] = -na*cf*alpha_f;
   left_eigenmatrix[0][2] = qs*bet2;
   left_eigenmatrix[0][3] = qs*bet3;
   left_eigenmatrix[0][4] = na*alpha_f/d;
   left_eigenmatrix[0][5] = as_prime*bet2;
   left_eigenmatrix[0][6] = as_prime*bet3;

/* left_eigenmatrix[1][0] = 0.0; */
/* left_eigenmatrix[1][1] = 0.0; */
   left_eigenmatrix[1][2] = -0.5*bet3;
   left_eigenmatrix[1][3] = 0.5*bet2;
/* left_eigenmatrix[1][4] = 0.0; */
   left_eigenmatrix[1][5] = -0.5*bet3*s/sqrtd;
   left_eigenmatrix[1][6] = 0.5*bet2*s/sqrtd;

/* left_eigenmatrix[2][0] = 0.0; */
   left_eigenmatrix[2][1] = -na*cs*alpha_s;
   left_eigenmatrix[2][2] = -qf*bet2;
   left_eigenmatrix[2][3] = -qf*bet3;
   left_eigenmatrix[2][4] = na*alpha_s/d;
   left_eigenmatrix[2][5] = -af_prime*bet2;
   left_eigenmatrix[2][6] = -af_prime*bet3;

   left_eigenmatrix[3][0] = 1.0;
/* left_eigenmatrix[3][1] = 0.0; */
/* left_eigenmatrix[3][2] = 0.0; */
/* left_eigenmatrix[3][3] = 0.0; */
   left_eigenmatrix[3][4] = -1.0/asq;
/* left_eigenmatrix[3][5] = 0.0; */
/* left_eigenmatrix[3][6] = 0.0; */

/* left_eigenmatrix[4][0] = 0.0; */
   left_eigenmatrix[4][1] = -left_eigenmatrix[2][1];
   left_eigenmatrix[4][2] = -left_eigenmatrix[2][2];
   left_eigenmatrix[4][3] = -left_eigenmatrix[2][3];
   left_eigenmatrix[4][4] = left_eigenmatrix[2][4];
   left_eigenmatrix[4][5] = left_eigenmatrix[2][5];
   left_eigenmatrix[4][6] = left_eigenmatrix[2][6];

/* left_eigenmatrix[5][0] = 0.0; */
/* left_eigenmatrix[5][1] = 0.0; */
   left_eigenmatrix[5][2] = -left_eigenmatrix[1][2];
   left_eigenmatrix[5][3] = -left_eigenmatrix[1][3];
/* left_eigenmatrix[5][4] = 0.0; */
   left_eigenmatrix[5][5] = left_eigenmatrix[1][5];
   left_eigenmatrix[5][6] = left_eigenmatrix[1][6];

/* left_eigenmatrix[6][0] = 0.0; */
   left_eigenmatrix[6][1] = -left_eigenmatrix[0][1];
   left_eigenmatrix[6][2] = -left_eigenmatrix[0][2];
   left_eigenmatrix[6][3] = -left_eigenmatrix[0][3];
   left_eigenmatrix[6][4] = left_eigenmatrix[0][4];
   left_eigenmatrix[6][5] = left_eigenmatrix[0][5];
   left_eigenmatrix[6][6] = left_eigenmatrix[0][6];
}

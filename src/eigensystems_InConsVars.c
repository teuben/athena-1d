/*==============================================================================
 *          File eigensystems_InConsVars (Multiple Functions)
 * This file contains functions to evaluate the eigenvalues, and left- and
 * right-eigenvectors of "Roe's matrix A" for the linearized system in the 
 * CONSERVED variables, i.e. U,t = AU,x, where U=(d,d*vx,d*vy,d*vz,[E],[By,Bz]).
 * The eigenvalues are returned through the argument list as a vector of length
 * NVAR.  The eigenvectors are returned as matrices of size (NVAR)x(NVAR), with
 * right-eigenvectors stored as COLUMNS (so R_i = right_eigenmatrix[*][i]), and
 * left-eigenvectors stored as ROWS (so L_i = left_eigenmatrix[i][*]).
 *
 * To improve performance components of the eigenvectors which are zero
 * are not set here (eigenmatrices must be initialized to zero in calling
 * routine).  However, for completeness, statements which set these values
 * are included but are commented out.
 *
 *============================================================================*/
#include <math.h>
#include <float.h>
#include "athena.def"
#include "athena.h"

#if defined(ISOTHERMAL) && defined(HYDRO)
/*---------------------------- ISOTHERMAL HYDRO --------------------------------
 * Input Arguments:
 *   v1,v2,v3 = Roe averaged components of velocity
 * Output Arguments:
 *   eigenvalues[NVAR],right_eigenmatrix[NVAR,NVAR],left_eigenmatrix[NVAR,NVAR]
 */
void Eigensystem4_isothermal_hydro_InConsVars(Real v1, Real v2, Real v3,
   Real eigenvalues[NVAR],
   Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR])
{

/* Compute eigenvalues (eq. B5) */

   eigenvalues[0] = v1 - ISOTHERMAL_C;
   eigenvalues[1] = v1;
   eigenvalues[2] = v1;
   eigenvalues[3] = v1 + ISOTHERMAL_C;

/* Right-eigenvectors, stored as COLUMNS (eq B3 sans 5th row and 4th column) */

   right_eigenmatrix[0][0] = 1.0;
   right_eigenmatrix[1][0] = v1 - ISOTHERMAL_C;
   right_eigenmatrix[2][0] = v2;
   right_eigenmatrix[3][0] = v3;

/* right_eigenmatrix[0][1] = 0.0; */
/* right_eigenmatrix[1][1] = 0.0; */
   right_eigenmatrix[2][1] = 1.0;
/* right_eigenmatrix[3][1] = 0.0; */

/* right_eigenmatrix[0][2] = 0.0; */
/* right_eigenmatrix[1][2] = 0.0; */
/* right_eigenmatrix[2][2] = 0.0; */
   right_eigenmatrix[3][2] = 1.0;

   right_eigenmatrix[0][3] = 1.0;
   right_eigenmatrix[1][3] = v1 + ISOTHERMAL_C;
   right_eigenmatrix[2][3] = v2;
   right_eigenmatrix[3][3] = v3;

/* Left-eigenvectors, stored as ROWS (eq. B7) */

   left_eigenmatrix[0][0] = 0.5*(1.0 + v1/ISOTHERMAL_C);
   left_eigenmatrix[0][1] = -0.5/ISOTHERMAL_C;
/* left_eigenmatrix[0][2] = 0.0; */
/* left_eigenmatrix[0][3] = 0.0; */

   left_eigenmatrix[1][0] = -v2;
/* left_eigenmatrix[1][1] = 0.0; */
   left_eigenmatrix[1][2] = 1.0;
/* left_eigenmatrix[1][3] = 0.0; */

   left_eigenmatrix[2][0] = -v3;
/* left_eigenmatrix[2][1] = 0.0; */
/* left_eigenmatrix[2][2] = 0.0; */
   left_eigenmatrix[2][3] = 1.0;

   left_eigenmatrix[3][0] = 0.5*(1.0 - v1/ISOTHERMAL_C);
   left_eigenmatrix[3][1] = 0.5/ISOTHERMAL_C;
/* left_eigenmatrix[3][2] = 0.0; */
/* left_eigenmatrix[3][3] = 0.0; */
}
#endif

#if defined(ADIABATIC) && defined(HYDRO)
/*------------------------------- ADIABATIC HYDRO ------------------------------
 *
 * Input: v1,v2,v3,h = Roe averaged velocities and enthalpy
 * Output: eigenvalues[NVAR]
 *         right_eigenmatrix[NVAR,NVAR], left_eigenmatrix[NVAR,NVAR];
 */
void Eigensystem4_adiabatic_hydro_InConsVars(Real v1, Real v2, Real v3, Real h,
   Real eigenvalues[NVAR],
   Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR])
{
  Real vsq,asq,a,na,qa;
  vsq = v1*v1 + v2*v2 + v3*v3;
  asq = GAMM1*MAX((h-0.5*vsq), TINY_NUMBER);
  a = sqrt(asq);

/* Compute eigenvalues (eq. B2) */

   eigenvalues[0] = v1 - a;
   eigenvalues[1] = v1;
   eigenvalues[2] = v1;
   eigenvalues[3] = v1;
   eigenvalues[4] = v1 + a;

/* Right-eigenvectors, stored as COLUMNS (eq. B3) */

   right_eigenmatrix[0][0] = 1.0;
   right_eigenmatrix[1][0] = v1 - a;
   right_eigenmatrix[2][0] = v2;
   right_eigenmatrix[3][0] = v3;
   right_eigenmatrix[4][0] = h - v1*a;

/* right_eigenmatrix[9][1] = 0.0; */
/* right_eigenmatrix[1][1] = 0.0; */
   right_eigenmatrix[2][1] = 1.0;
/* right_eigenmatrix[3][1] = 0.0; */
   right_eigenmatrix[4][1] = v2;

/* right_eigenmatrix[0][2] = 0.0; */
/* right_eigenmatrix[1][2] = 0.0; */
/* right_eigenmatrix[2][2] = 0.0; */
   right_eigenmatrix[3][2] = 1.0;
   right_eigenmatrix[4][2] = v3;

   right_eigenmatrix[0][3] = 1.0;
   right_eigenmatrix[1][3] = v1;
   right_eigenmatrix[2][3] = v2;
   right_eigenmatrix[3][3] = v3;
   right_eigenmatrix[4][3] = 0.5*vsq;

   right_eigenmatrix[0][4] = 1.0;
   right_eigenmatrix[1][4] = v1 + a;
   right_eigenmatrix[2][4] = v2;
   right_eigenmatrix[3][4] = v3;
   right_eigenmatrix[4][4] = h + v1*a;

/* Left-eigenvectors, stored as ROWS (eq. B4) */

   na = 0.5/asq;
   left_eigenmatrix[0][0] = na*(0.5*GAMM1*vsq + v1*a);
   left_eigenmatrix[0][1] = -na*(GAMM1*v1 + a);
   left_eigenmatrix[0][2] = -na*GAMM1*v2;
   left_eigenmatrix[0][3] = -na*GAMM1*v3;
   left_eigenmatrix[0][4] = na*GAMM1;

   left_eigenmatrix[1][0] = -v2;
/* left_eigenmatrix[1][1] = 0.0; */
   left_eigenmatrix[1][2] = 1.0;
/* left_eigenmatrix[1][3] = 0.0; */
/* left_eigenmatrix[1][4] = 0.0; */

   left_eigenmatrix[2][0] = -v3;
/* left_eigenmatrix[2][1] = 0.0; */
/* left_eigenmatrix[2][2] = 0.0; */
   left_eigenmatrix[2][3] = 1.0;
/* left_eigenmatrix[2][4] = 0.0; */

   qa = GAMM1/asq;
   left_eigenmatrix[3][0] = 1.0 - na*GAMM1*vsq;
   left_eigenmatrix[3][1] = qa*v1;
   left_eigenmatrix[3][2] = qa*v2;
   left_eigenmatrix[3][3] = qa*v3;
   left_eigenmatrix[3][4] = -qa;

   left_eigenmatrix[4][0] = na*(0.5*GAMM1*vsq - v1*a);
   left_eigenmatrix[4][1] = -na*(GAMM1*v1 - a);
   left_eigenmatrix[4][2] = left_eigenmatrix[0][2];
   left_eigenmatrix[4][3] = left_eigenmatrix[0][3];
   left_eigenmatrix[4][4] = left_eigenmatrix[0][4];
}
#endif


#if defined(ISOTHERMAL) && defined(MHD)
/*------------------------------ ISOTHERMAL MHD --------------------------------
 *
 * Input: d,v1,v2,v3,b1,b2,b3 = Roe averaged density, velocities, and B field
 *        x,y = numerical factors (eqs B15 and B16)
 * Output: eigenvalues[NVAR]
 *         right_eigenmatrix[NVAR,NVAR], left_eigenmatrix[NVAR,NVAR];
 */
void Eigensystem4_isothermal_mhd_InConsVars(Real d, Real v1, Real v2, Real v3,
   Real b1, Real b2, Real b3, Real x, Real y, Real eigenvalues[NVAR],
   Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR])
{
  Real btsq,bt_starsq,vaxsq,twid_csq,q_starsq,cfsq,cf,cssq,cs;
  Real bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,alpha_f,alpha_s;
  Real sqrtd,s,twid_c,qf,qs,af_prime,as_prime,vax;
  Real norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
  btsq = b2*b2 + b3*b3;
  bt_starsq = btsq*y;

/* Compute fast- and slow-magnetosonic speeds (eq. B42 and B43) */

   vaxsq = b1*b1/d;
   twid_csq = (ISOTHERMAL_C_SQ) + x;
   q_starsq  = twid_csq + (vaxsq + bt_starsq/d);
   cfsq = 0.5*(q_starsq + sqrt(q_starsq*q_starsq - 4.0*twid_csq*vaxsq));
   cf  = sqrt(cfsq);
   cssq = 0.5*(q_starsq - sqrt(q_starsq*q_starsq - 4.0*twid_csq*vaxsq));
   if (cssq <= 0.0) {
      cssq = 0.0;
   }
   cs  = sqrt(cssq);

/* Compute beta's (eqs A13, B29) */

   bt = sqrt(btsq);
   bt_star = sqrt(bt_starsq);
   if (bt == 0.0) {
      bet2 = ONE_OVER_SQRT2 ;
      bet3 = ONE_OVER_SQRT2 ;
      bet2_star = ONE_OVER_SQRT2 ;
      bet3_star = ONE_OVER_SQRT2 ;
   } 
   else {
      bet2 = b2/bt;
      bet3 = b3/bt;
      bet2_star = b2/bt_star;
      bet3_star = b3/bt_star;
   }
   bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;

/* Compute alpha's (eq A12) */

   if ((cfsq-cssq) == 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
   } else if ((twid_csq - cssq) <= 0.0) {
      alpha_f = 0.0;
      alpha_s = 1.0;
   } else if ((cfsq - twid_csq) <= 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
   } else {
      alpha_f = sqrt((twid_csq - cssq)/(cfsq - cssq));
      alpha_s = sqrt((cfsq - twid_csq)/(cfsq - cssq));
   }

/* Compute Q's (eq. A11), etc. */

   sqrtd = sqrt(d);
   s = SIGN(b1);
   twid_c = sqrt(twid_csq);
   qf = cf*alpha_f*s;
   qs = cs*alpha_s*s;
   af_prime = twid_c*alpha_f/sqrtd;
   as_prime = twid_c*alpha_s/sqrtd;

/* Compute eigenvalues (eq. XX) */

   vax  = sqrt(vaxsq);
   eigenvalues[0] = v1 - cf;
   eigenvalues[1] = v1 - vax;
   eigenvalues[2] = v1 - cs;
   eigenvalues[3] = v1 + cs;
   eigenvalues[4] = v1 + vax;
   eigenvalues[5] = v1 + cf;

/* Right-eigenvectors, stored as COLUMNS (eq B22 sans 5th row and 4th column) */

   right_eigenmatrix[0][0] = alpha_f;
   right_eigenmatrix[1][0] = alpha_f*(v1 - cf);
   right_eigenmatrix[2][0] = alpha_f*v2 + qs*bet2_star;
   right_eigenmatrix[3][0] = alpha_f*v3 + qs*bet3_star;
   right_eigenmatrix[4][0] = as_prime*bet2_star;
   right_eigenmatrix[5][0] = as_prime*bet3_star;

/* right_eigenmatrix[0][1] = 0.0; */
/* right_eigenmatrix[1][1] = 0.0; */
   right_eigenmatrix[2][1] = -bet3;
   right_eigenmatrix[3][1] = bet2;
   right_eigenmatrix[4][1] = -bet3*s/sqrtd;
   right_eigenmatrix[5][1] = bet2*s/sqrtd;

   right_eigenmatrix[0][2] = alpha_s;
   right_eigenmatrix[1][2] = alpha_s*(v1 - cs);
   right_eigenmatrix[2][2] = alpha_s*v2 - qf*bet2_star;
   right_eigenmatrix[3][2] = alpha_s*v3 - qf*bet3_star;
   right_eigenmatrix[4][2] = -af_prime*bet2_star;
   right_eigenmatrix[5][2] = -af_prime*bet3_star;

   right_eigenmatrix[0][3] = alpha_s;
   right_eigenmatrix[1][3] = alpha_s*(v1 + cs);
   right_eigenmatrix[2][3] = alpha_s*v2 + qf*bet2_star;
   right_eigenmatrix[3][3] = alpha_s*v3 + qf*bet3_star;
   right_eigenmatrix[4][3] = right_eigenmatrix[4][2];
   right_eigenmatrix[5][3] = right_eigenmatrix[5][2];

/* right_eigenmatrix[0][4] = 0.0; */
/* right_eigenmatrix[1][4] = 0.0; */
   right_eigenmatrix[2][4] = bet3;
   right_eigenmatrix[3][4] = -bet2;
   right_eigenmatrix[4][4] = right_eigenmatrix[4][1];
   right_eigenmatrix[5][4] = right_eigenmatrix[5][1];

   right_eigenmatrix[0][5] = alpha_f;
   right_eigenmatrix[1][5] = alpha_f*(v1 + cf);
   right_eigenmatrix[2][5] = alpha_f*v2 - qs*bet2_star;
   right_eigenmatrix[3][5] = alpha_f*v3 - qs*bet3_star;
   right_eigenmatrix[4][5] = right_eigenmatrix[4][0];
   right_eigenmatrix[5][5] = right_eigenmatrix[5][0];

/* Left-eigenvectors, stored as ROWS (eq. B44) */

/* Normalize by 1/2a^{2}: quantities denoted by \hat{f} */
   norm = 0.5/twid_csq;
   cff = norm*alpha_f*cf;
   css = norm*alpha_s*cs;
   qf *= norm;
   qs *= norm;
   af = norm*af_prime*d;
   as = norm*as_prime*d;
   afpb = norm*af_prime*bt_star;
   aspb = norm*as_prime*bt_star;

   q2_star = bet2_star/bet_starsq;
   q3_star = bet3_star/bet_starsq;
   vqstr = (v2*q2_star + v3*q3_star);

   left_eigenmatrix[0][0] = cff*(cf+v1) - qs*vqstr - aspb;
   left_eigenmatrix[0][1] = -cff;
   left_eigenmatrix[0][2] = qs*q2_star;
   left_eigenmatrix[0][3] = qs*q3_star;
   left_eigenmatrix[0][4] = as*q2_star;
   left_eigenmatrix[0][5] = as*q3_star;

   left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
/* left_eigenmatrix[1][1] = 0.0; */
   left_eigenmatrix[1][2] = -0.5*bet3;
   left_eigenmatrix[1][3] = 0.5*bet2;
   left_eigenmatrix[1][4] = -0.5*sqrtd*bet3*s;
   left_eigenmatrix[1][5] = 0.5*sqrtd*bet2*s;

   left_eigenmatrix[2][0] = css*(cs+v1) + qf*vqstr + afpb;
   left_eigenmatrix[2][1] = -css;
   left_eigenmatrix[2][2] = -qf*q2_star;
   left_eigenmatrix[2][3] = -qf*q3_star;
   left_eigenmatrix[2][4] = -af*q2_star;
   left_eigenmatrix[2][5] = -af*q3_star;

   left_eigenmatrix[3][0] = css*(cs-v1) - qf*vqstr + afpb;
   left_eigenmatrix[3][1] = css;
   left_eigenmatrix[3][2] = -left_eigenmatrix[2][2];
   left_eigenmatrix[3][3] = -left_eigenmatrix[2][3];
   left_eigenmatrix[3][4] = left_eigenmatrix[2][4];
   left_eigenmatrix[3][5] = left_eigenmatrix[2][5];

   left_eigenmatrix[4][0] = -left_eigenmatrix[1][0];
/* left_eigenmatrix[4][1] = 0.0; */
   left_eigenmatrix[4][2] = -left_eigenmatrix[1][2];
   left_eigenmatrix[4][3] = -left_eigenmatrix[1][3];
   left_eigenmatrix[4][4] = left_eigenmatrix[1][4];
   left_eigenmatrix[4][5] = left_eigenmatrix[1][5];

   left_eigenmatrix[5][0] = cff*(cf-v1) + qs*vqstr - aspb;
   left_eigenmatrix[5][1] = cff;
   left_eigenmatrix[5][2] = -left_eigenmatrix[0][2];
   left_eigenmatrix[5][3] = -left_eigenmatrix[0][3];
   left_eigenmatrix[5][4] = left_eigenmatrix[0][4];
   left_eigenmatrix[5][5] = left_eigenmatrix[0][5];
}
#endif

#if defined(ADIABATIC) && defined(MHD)
/*------------------------------- ADIABATIC MHD --------------------------------
 *
 * Input: d,v1,v2,v3,h,b1,b2,b3 = Roe averaged density, velocities, enthalpy, B
 *        x,y = numerical factors (see eqn XX)
 * Output: eigenvalues[NVAR]
 *         right_eigenmatrix[NVAR,NVAR], left_eigenmatrix[NVAR,NVAR];
 */
void Eigensystem4_adiabatic_mhd_InConsVars(Real d, Real v1, Real v2, Real v3,
   Real h, Real b1, Real b2, Real b3, Real x, Real y,
   Real eigenvalues[NVAR],
   Real right_eigenmatrix[NVAR][NVAR], Real left_eigenmatrix[NVAR][NVAR])
{
  Real vsq,btsq,bt_starsq,vaxsq,hp,twid_asq,q_starsq,cfsq,cf,cssq,cs;
  Real bt,bt_star,bet2,bet3,bet2_star,bet3_star,bet_starsq,vbet,alpha_f,alpha_s;
  Real sqrtd,s,twid_a,qf,qs,af_prime,as_prime,afpbb,aspbb,vax;
  Real norm,cff,css,af,as,afpb,aspb,q2_star,q3_star,vqstr;
  vsq = v1*v1 + v2*v2 + v3*v3;
  btsq = b2*b2 + b3*b3;
  bt_starsq = (GAMM1 - GAMM2*y)*btsq;

/* Compute fast- and slow-magnetosonic speeds (eqs. B18-20) */

   vaxsq = b1*b1/d;
   hp = h - (vaxsq + btsq/d);
   twid_asq = MAX((GAMM1*(hp-0.5*vsq)-GAMM2*x), TINY_NUMBER);
   q_starsq  = twid_asq + (vaxsq + bt_starsq/d);
   cfsq = 0.5*(q_starsq + sqrt(q_starsq*q_starsq - 4.0*twid_asq*vaxsq));
   cf  = sqrt(cfsq);
   cssq = 0.5*(q_starsq - sqrt(q_starsq*q_starsq - 4.0*twid_asq*vaxsq));
   if (cssq <= 0.0) {
      cssq = 0.0;
   }
   cs  = sqrt(cssq);

/* Compute beta(s) (eqs A13, B28) */

   bt = sqrt(btsq);
   bt_star = sqrt(bt_starsq);
   if (bt == 0.0) {
      bet2 = ONE_OVER_SQRT2 ;
      bet3 = ONE_OVER_SQRT2 ;
      bet2_star = ONE_OVER_SQRT2 ;
      bet3_star = ONE_OVER_SQRT2 ;
   } else {
      bet2 = b2/bt;
      bet3 = b3/bt;
      bet2_star = b2/bt_star;
      bet3_star = b3/bt_star;
   }
   bet_starsq = bet2_star*bet2_star + bet3_star*bet3_star;
   vbet = v2*bet2_star + v3*bet3_star;

/* Compute alpha(s) (eq A12) */

   if ((cfsq-cssq) == 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
   } else if ( (twid_asq - cssq) <= 0.0) {
      alpha_f = 0.0;
      alpha_s = 1.0;
   } else if ( (cfsq - twid_asq) <= 0.0) {
      alpha_f = 1.0;
      alpha_s = 0.0;
   } else {
      alpha_f = sqrt((twid_asq - cssq)/(cfsq - cssq));
      alpha_s = sqrt((cfsq - twid_asq)/(cfsq - cssq));
   }

/* Compute Q(s) and A(s) (eq. A11), etc. */

   sqrtd = sqrt(d);
   s = SIGN(b1);
   twid_a = sqrt(twid_asq);
   qf = cf*alpha_f*s;
   qs = cs*alpha_s*s;
   af_prime = twid_a*alpha_f/sqrtd;
   as_prime = twid_a*alpha_s/sqrtd;
   afpbb = af_prime*bt_star*bet_starsq;
   aspbb = as_prime*bt_star*bet_starsq;

/* Compute eigenvalues (eq. B17) */

   vax = sqrt(vaxsq);
   eigenvalues[0] = v1 - cf;
   eigenvalues[1] = v1 - vax;
   eigenvalues[2] = v1 - cs;
   eigenvalues[3] = v1;
   eigenvalues[4] = v1 + cs;
   eigenvalues[5] = v1 + vax;
   eigenvalues[6] = v1 + cf;

/* Right-eigenvectors, stored as COLUMNS (eq B21) */

   right_eigenmatrix[0][0] = alpha_f;
   right_eigenmatrix[1][0] = alpha_f*(v1 - cf);
   right_eigenmatrix[2][0] = alpha_f*v2 + qs*bet2_star;
   right_eigenmatrix[3][0] = alpha_f*v3 + qs*bet3_star;
   right_eigenmatrix[4][0] = alpha_f*(hp - v1*cf) + qs*vbet + aspbb;
   right_eigenmatrix[5][0] = as_prime*bet2_star;
   right_eigenmatrix[6][0] = as_prime*bet3_star;

/* right_eigenmatrix[0][1] = 0.0; */
/* right_eigenmatrix[1][1] = 0.0; */
   right_eigenmatrix[2][1] = -bet3;
   right_eigenmatrix[3][1] = bet2;
   right_eigenmatrix[4][1] = -(v2*bet3 - v3*bet2);
   right_eigenmatrix[5][1] = -bet3*s/sqrtd;
   right_eigenmatrix[6][1] = bet2*s/sqrtd;

   right_eigenmatrix[0][2] = alpha_s;
   right_eigenmatrix[1][2] = alpha_s*(v1 - cs);
   right_eigenmatrix[2][2] = alpha_s*v2 - qf*bet2_star;
   right_eigenmatrix[3][2] = alpha_s*v3 - qf*bet3_star;
   right_eigenmatrix[4][2] = alpha_s*(hp - v1*cs) - qf*vbet - afpbb;
   right_eigenmatrix[5][2] = -af_prime*bet2_star;
   right_eigenmatrix[6][2] = -af_prime*bet3_star;

   right_eigenmatrix[0][3] = 1.0;
   right_eigenmatrix[1][3] = v1;
   right_eigenmatrix[2][3] = v2;
   right_eigenmatrix[3][3] = v3;
   right_eigenmatrix[4][3] = 0.5*vsq + GAMM2*x/GAMM1;
/* right_eigenmatrix[5][3] = 0.0; */
/* right_eigenmatrix[6][3] = 0.0; */

   right_eigenmatrix[0][4] = alpha_s;
   right_eigenmatrix[1][4] = alpha_s*(v1 + cs);
   right_eigenmatrix[2][4] = alpha_s*v2 + qf*bet2_star;
   right_eigenmatrix[3][4] = alpha_s*v3 + qf*bet3_star;
   right_eigenmatrix[4][4] = alpha_s*(hp + v1*cs) + qf*vbet - afpbb;
   right_eigenmatrix[5][4] = right_eigenmatrix[5][2];
   right_eigenmatrix[6][4] = right_eigenmatrix[6][2];

/* right_eigenmatrix[0][5] = 0.0; */
/* right_eigenmatrix[1][5] = 0.0; */
   right_eigenmatrix[2][5] = bet3;
   right_eigenmatrix[3][5] = -bet2;
   right_eigenmatrix[4][5] = -right_eigenmatrix[4][1];
   right_eigenmatrix[5][5] = right_eigenmatrix[5][1];
   right_eigenmatrix[6][5] = right_eigenmatrix[6][1];

   right_eigenmatrix[0][6] = alpha_f;
   right_eigenmatrix[1][6] = alpha_f*(v1 + cf);
   right_eigenmatrix[2][6] = alpha_f*v2 - qs*bet2_star;
   right_eigenmatrix[3][6] = alpha_f*v3 - qs*bet3_star;
   right_eigenmatrix[4][6] = alpha_f*(hp + v1*cf) - qs*vbet + aspbb;
   right_eigenmatrix[5][6] = right_eigenmatrix[5][0];
   right_eigenmatrix[6][6] = right_eigenmatrix[6][0];

/* Left-eigenvectors, stored as ROWS (eq. B30) */

/* Normalize by 1/2a^{2}: quantities denoted by \hat{f} */
   norm = 0.5/twid_asq;
   cff = norm*alpha_f*cf;
   css = norm*alpha_s*cs;
   qf *= norm;
   qs *= norm;
   af = norm*af_prime*d;
   as = norm*as_prime*d;
   afpb = norm*af_prime*bt_star;
   aspb = norm*as_prime*bt_star;

/* Normalize by (gamma-1)/2a^{2}: quantities denoted by \bar{f} */
   norm *= GAMM1;
   alpha_f *= norm;
   alpha_s *= norm;
   q2_star = bet2_star/bet_starsq;
   q3_star = bet3_star/bet_starsq;
   vqstr = (v2*q2_star + v3*q3_star);
   norm *= 2.0;

   left_eigenmatrix[0][0] = alpha_f*(vsq-hp) + cff*(cf+v1) - qs*vqstr - aspb;
   left_eigenmatrix[0][1] = -alpha_f*v1 - cff;
   left_eigenmatrix[0][2] = -alpha_f*v2 + qs*q2_star;
   left_eigenmatrix[0][3] = -alpha_f*v3 + qs*q3_star;
   left_eigenmatrix[0][4] = alpha_f;
   left_eigenmatrix[0][5] = as*q2_star - alpha_f*b2;
   left_eigenmatrix[0][6] = as*q3_star - alpha_f*b3;

   left_eigenmatrix[1][0] = 0.5*(v2*bet3 - v3*bet2);
/* left_eigenmatrix[1][1] = 0.0; */
   left_eigenmatrix[1][2] = -0.5*bet3;
   left_eigenmatrix[1][3] = 0.5*bet2;
/* left_eigenmatrix[1][4] = 0.0; */
   left_eigenmatrix[1][5] = -0.5*sqrtd*bet3*s;
   left_eigenmatrix[1][6] = 0.5*sqrtd*bet2*s;

   left_eigenmatrix[2][0] = alpha_s*(vsq-hp) + css*(cs+v1) + qf*vqstr + afpb;
   left_eigenmatrix[2][1] = -alpha_s*v1 - css;
   left_eigenmatrix[2][2] = -alpha_s*v2 - qf*q2_star;
   left_eigenmatrix[2][3] = -alpha_s*v3 - qf*q3_star;
   left_eigenmatrix[2][4] = alpha_s;
   left_eigenmatrix[2][5] = -af*q2_star - alpha_s*b2;
   left_eigenmatrix[2][6] = -af*q3_star - alpha_s*b3;

   left_eigenmatrix[3][0] = 1.0 - norm*(0.5*vsq + GAMM2*x/GAMM1); 
   left_eigenmatrix[3][1] = norm*v1;
   left_eigenmatrix[3][2] = norm*v2;
   left_eigenmatrix[3][3] = norm*v3;
   left_eigenmatrix[3][4] = -norm;
   left_eigenmatrix[3][5] = norm*b2;
   left_eigenmatrix[3][6] = norm*b3;

   left_eigenmatrix[4][0] = alpha_s*(vsq-hp) + css*(cs-v1) - qf*vqstr + afpb;
   left_eigenmatrix[4][1] = -alpha_s*v1 + css;
   left_eigenmatrix[4][2] = -alpha_s*v2 + qf*q2_star;
   left_eigenmatrix[4][3] = -alpha_s*v3 + qf*q3_star;
   left_eigenmatrix[4][4] = alpha_s;
   left_eigenmatrix[4][5] = left_eigenmatrix[2][5];
   left_eigenmatrix[4][6] = left_eigenmatrix[2][6];

   left_eigenmatrix[5][0] = -left_eigenmatrix[1][0];
/* left_eigenmatrix[5][1] = 0.0; */
   left_eigenmatrix[5][2] = -left_eigenmatrix[1][2];
   left_eigenmatrix[5][3] = -left_eigenmatrix[1][3];
/* left_eigenmatrix[5][4] = 0.0; */
   left_eigenmatrix[5][5] = left_eigenmatrix[1][5];
   left_eigenmatrix[5][6] = left_eigenmatrix[1][6];

   left_eigenmatrix[6][0] = alpha_f*(vsq-hp) + cff*(cf-v1) + qs*vqstr - aspb;
   left_eigenmatrix[6][1] = -alpha_f*v1 + cff;
   left_eigenmatrix[6][2] = -alpha_f*v2 - qs*q2_star;
   left_eigenmatrix[6][3] = -alpha_f*v3 - qs*q3_star;
   left_eigenmatrix[6][4] = alpha_f;
   left_eigenmatrix[6][5] = left_eigenmatrix[0][5];
   left_eigenmatrix[6][6] = left_eigenmatrix[0][6];
}
#endif

/*==============================================================================
 *                           Function SET_BVAL_ARRAYS
 *
 *============================================================================*/
#include <stdio.h>
#include "athena.def"
#include "athena.h"
void set_bval_arrays(struct grid_block *agrid, struct bval_array *bval)
{
int i=0,n,ii,is,ie;
REAL reflect_s;
   is = agrid->is; ie = agrid->ie;

/*------------------------ Start with Inner I boundary -----------------------*/

   if (agrid->niib == 1) {reflect_s = -1.0;} else {reflect_s = 1.0;}

/* Do niib=3 (inflow) */

   if (agrid->niib == 3) {
      for (i=0;  i<=3;  i++) {
      for (n=0;  n<NVAR;  n++) {
         bval->uiib[i][n] = agrid->boundary_values.uiib[i][n];
         bval->bxiib[i] = agrid->boundary_values.bxiib[i];
      }}

/* Do niib=1 (reflection), 2 (flow-out), 4 (periodic) by adjusting i-index   */

   } else {
      for (i=0;  i<=3;  i++) {
         if (agrid->niib == 1) {ii = is+(3-i);}
         if (agrid->niib == 2) {ii = is;}
         if (agrid->niib == 4) {ii = ie-(3-i);}
         for (n=0;  n<NVAR;  n++) {
            bval->uiib[i][n] = agrid->u[ii][n];
            bval->bxiib[i] = agrid->bx[ii];
         }
         bval->uiib[i][1] = agrid->u[ii][1]*reflect_s;
      }
   }

/*----------------------- Continue with Outer I boundary ---------------------*/

   if (agrid->noib == 1) {reflect_s = -1.0;} else {reflect_s = 1.0;}

/* Do noib=3 (inflow) */

   if (agrid->noib == 3) {
      for (i=0;  i<=3;  i++) {
      for (n=0;  n<NVAR;  n++) {
         bval->uoib[i][n] = agrid->boundary_values.uoib[i][n];
         bval->bxoib[i] = agrid->boundary_values.bxoib[i];
      }}

/* Do noib=1 (reflection), 2 (flow-out), 4 (periodic) by adjusting i-index   */
/* Note Bx ghost zones start at ie+2 rather than ie+1 */

   } else {
      for (i=0;  i<=3;  i++) {
         if (agrid->noib == 1) {ii = ie-i;}
         if (agrid->noib == 2) {ii = ie;}
         if (agrid->noib == 4) {ii = is+i;}
         for (n=0;  n<NVAR;  n++) {
            bval->uoib[i][n] = agrid->u[ii][n];
            bval->bxoib[i] = agrid->bx[ii+1];
         }
         bval->uoib[i][1] = agrid->u[ii][1]*reflect_s;
      }
   }
}

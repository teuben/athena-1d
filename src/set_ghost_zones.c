/*==============================================================================
 *                        Function SET_GHOST_ZONES 
 *
 *============================================================================*/
#include "athena.def"
#include "athena.h"

void set_ghost_zones(struct grid_block *agrid, struct bval_array *bval)
{
  int i=0,n;
  int is,ie;

  is = agrid->is; ie = agrid->ie;

/*========================= Inner/Outer I boundaries =========================*/

  for (i=0;  i<=3;   i++) {
    for (n=0;  n<NVAR; n++) {
      agrid->u[i][n] = bval->uiib[i][n];
      agrid->u[ie+(i+1)][n] = bval->uoib[i][n];
    }
    agrid->bx[i] = bval->bxiib[i];
    agrid->bx[ie+(i+1)] = bval->bxoib[i];
  }
}

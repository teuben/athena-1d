/*==============================================================================
 *                        Function SET_GHOST_ZONES 
 *
 *============================================================================*/
#include "athena.def"
#include "athena.h"

#if defined ONE_D

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

#elif defined TWO_D

void set_ghost_zones(struct grid_block *agrid, struct bval_array *bval){
  int i,j,n;
  int is1,is2,ie1,ie2;

  is1 = agrid->is1;  ie1 = agrid->ie1;
  is2 = agrid->is2;  ie2 = agrid->ie2;

/*========================= Inner/Outer 1 boundaries =========================*/

  for(i=0; i<4; i++){
    for(j=is2-4; j<=ie2+4; j++){
      for(n=0; n<NVAR; n++){
	agrid->u[i][j][n] = bval->uiib1[i][j][n];
	agrid->u[ie1+(i+1)][j][n] = bval->uoib1[i][j][n];
      }
      agrid->bx[i][j] = bval->bxiib1[i][j];
      agrid->bx[ie1+(i+1)][j] = bval->bxoib1[i][j];

      agrid->by[i][j] = bval->byiib1[i][j];
      agrid->by[ie1+(i+1)][j] = bval->byoib1[i][j];
    }
  }

/*========================= Inner/Outer 2 boundaries =========================*/

  for(i=is1-4; i<=ie1+4; i++){
    for(j=0; j<4; j++){
      for(n=0; n<NVAR; n++){
	agrid->u[i][j][n] = bval->uiib2[i][j][n];
	agrid->u[i][ie2+(j+1)][n] = bval->uoib2[i][j][n];
      }
      agrid->bx[i][j] = bval->bxiib2[i][j];
      agrid->bx[i][ie2+(j+1)] = bval->bxoib2[i][j];

      agrid->by[i][j] = bval->byiib2[i][j];
      agrid->by[i][ie2+(j+1)] = bval->byoib2[i][j];
    }
  }
}

#else /* THREE_D */


#endif

/*==============================================================================
 *                           Function SET_BVAL_ARRAYS
 *
 *============================================================================*/
#include <stdio.h>
#include "athena.def"
#include "athena.h"


#if defined ONE_D


void set_bval_arrays(struct grid_block *agrid, struct bval_array *bval)
{
  int i=0,n,ii,is,ie;
  Real reflect_s;

  is = agrid->is; ie = agrid->ie;

/*------------------------ Start with Inner I boundary -----------------------*/

  if (agrid->niib == 1) {reflect_s = -1.0;} else {reflect_s = 1.0;}

/* Do niib=3 (inflow) */

  if (agrid->niib == 3) {
    for (i=0;  i<=3;  i++) {
      for (n=0;  n<NVAR;  n++) {
	bval->uiib[i][n] = agrid->boundary_values.uiib[i][n];
      }
#ifdef MHD
      bval->bxiib[i] = agrid->boundary_values.bxiib[i];
#endif /* MHD */
    }

/* Do niib=1 (reflection), 2 (flow-out), 4 (periodic) by adjusting i-index   */

  } else {
    for (i=0;  i<=3;  i++) {
      switch(agrid->niib){
      case 1:
	ii = is+(3-i);  break;
      case 2:
	ii = is;  break;
      case 4:
	ii = ie-(3-i);  break;
      default:
	fprintf(stderr,"[set_bval_arrays]: agrid->niib = %d unknown\n",
		agrid->niib);
	ii = is; /* Choose outflow as default */
      }
      for (n=0;  n<NVAR;  n++) {
	bval->uiib[i][n] = agrid->u[ii][n];
      }
#ifdef MHD
      bval->bxiib[i] = agrid->bx[ii];
#endif /* MHD */
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
      }
#ifdef MHD
      bval->bxoib[i] = agrid->boundary_values.bxoib[i];
#endif /* MHD */
    }

/* Do noib=1 (reflection), 2 (flow-out), 4 (periodic) by adjusting i-index   */
/* Note Bx ghost zones start at ie+2 rather than ie+1 */

  } else {
    for (i=0;  i<=3;  i++) {
      switch(agrid->noib){
      case 1:
	ii = ie-i;  break;
      case 2:
	ii = ie;  break;
      case 4:
	ii = is+i;  break;
      default:
	fprintf(stderr,"[set_bval_arrays]: agrid->noib = %d unknown\n",
		agrid->noib);
	ii = ie; /* Choose outflow as default */
      }
      for (n=0;  n<NVAR;  n++) {
	bval->uoib[i][n] = agrid->u[ii][n];
      }
#ifdef MHD
      bval->bxoib[i] = agrid->bx[ii+1];
#endif /* MHD */
      bval->uoib[i][1] = agrid->u[ii][1]*reflect_s;
    }
  }
}


#elif defined TWO_D


void set_bval_arrays(struct grid_block *agrid, struct bval_array *bval){
  int i,j,n,ii,jj;
  int is1,is2,ie1,ie2;
  Real reflect_s;

  is1 = agrid->is1;  ie1 = agrid->ie1;
  is2 = agrid->is2;  ie2 = agrid->ie2;

/*------------------------ Start with Inner 1 boundary -----------------------*/

  if (agrid->niib1 == 1) {reflect_s = -1.0;} else {reflect_s = 1.0;}

/* Do niib1=3 (inflow) */

  if (agrid->niib1 == 3) {
    for (i=0; i<4; i++) {
      for (j=is2-4; j<=ie2+4; j++){
	for (n=0; n<NVAR; n++) {
	  bval->uiib1[i][j][n] = agrid->boundary_values.uiib1[i][j][n];
	}
#ifdef MHD
	bval->bxiib1[i][j] = agrid->boundary_values.bxiib1[i][j];
	bval->byiib1[i][j] = agrid->boundary_values.byiib1[i][j];
#endif /* MHD */
      }
    }

/* Do niib1=1 (reflection), 2 (flow-out), 4 (periodic) by adjusting i-index   */

  } else {
    for (i=0; i<4; i++) {
      switch(agrid->niib1){
      case 1:
	ii = is1+(3-i);  break;
      case 2:
	ii = is1;  break;
      case 4:
	ii = ie1-(3-i);  break;
      default:
	fprintf(stderr,"[set_bval_arrays]: agrid->niib1 = %d unknown\n",
		agrid->niib1);
	ii = is1; /* Choose outflow as default */
      }
      for (j=is2-4; j<=ie2+4; j++){
	for (n=0; n<NVAR; n++) {
	  bval->uiib1[i][j][n] = agrid->u[ii][j][n];
	}
#ifdef MHD
	bval->bxiib1[i][j] = agrid->bx[ii][j];
	bval->byiib1[i][j] = agrid->by[ii][j];
#endif /* MHD */
	bval->uiib1[i][j][1] = agrid->u[ii][j][1]*reflect_s;
      }
    }
  }

/*----------------------- Continue with Outer 1 boundary ---------------------*/

  if (agrid->noib1 == 1) {reflect_s = -1.0;} else {reflect_s = 1.0;}

/* Do noib1=3 (inflow) */

  if (agrid->noib1 == 3) {
    for (i=0; i<4; i++) {
      for (j=is2-4; j<=ie2+4; j++){
	for (n=0; n<NVAR; n++) {
	  bval->uoib1[i][j][n] = agrid->boundary_values.uoib1[i][j][n];
	}
#ifdef MHD
	bval->bxoib1[i][j] = agrid->boundary_values.bxoib1[i][j];
	bval->byoib1[i][j] = agrid->boundary_values.byoib1[i][j];
#endif /* MHD */
      }
    }

/* Do noib1=1 (reflection), 2 (flow-out), 4 (periodic) by adjusting i-index   */
/* Note Bx ghost zones start at ie+2 rather than ie+1 */

  } else {
    for (i=0; i<4; i++) {
      switch(agrid->noib1){
      case 1:
	ii = ie1-i;  break;
      case 2:
	ii = ie1;  break;
      case 4:
	ii = is1+i;  break;
      default:
	fprintf(stderr,"[set_bval_arrays]: agrid->noib1 = %d unknown\n",
		agrid->noib1);
	ii = ie1; /* Choose outflow as default */
      }
      for(j=is2-4; j<=ie2+4; j++){
	for (n=0; n<NVAR; n++) {
	  bval->uoib1[i][j][n] = agrid->u[ii][j][n];
	}
#ifdef MHD
	bval->bxoib1[i][j] = agrid->bx[ii+1][j];
	bval->byoib1[i][j] = agrid->by[ii][j];
#endif /* MHD */
	bval->uoib1[i][j][1] = agrid->u[ii][j][1]*reflect_s;
      }
    }
  }


/*----------------------- Continue with Inner 2 boundary ---------------------*/

  if (agrid->niib2 == 1) {reflect_s = -1.0;} else {reflect_s = 1.0;}

/* Do niib2=3 (inflow) */

  if (agrid->niib2 == 3) {
    for (i=is1-4; i<=ie1+4; i++){
      for (j=0; j<4; j++) {
	for (n=0; n<NVAR; n++) {
	  bval->uiib2[i][j][n] = agrid->boundary_values.uiib2[i][j][n];
	}
#ifdef MHD
	bval->bxiib2[i][j] = agrid->boundary_values.bxiib2[i][j];
	bval->byiib2[i][j] = agrid->boundary_values.byiib2[i][j];
#endif /* MHD */
      }
    }

/* Do niib2=1 (reflection), 2 (flow-out), 4 (periodic) by adjusting i-index   */

  } else {
    for (i=is1-4; i<=ie1+4; i++){
      for (j=0; j<4; j++) {
	switch(agrid->niib2){
	case 1:
	  jj = is2+(3-j);  break;
	case 2:
	  jj = is2;  break;
	case 4:
	  jj = ie2-(3-j);  break;
	default:
	  fprintf(stderr,"[set_bval_arrays]: agrid->niib2 = %d unknown\n",
		  agrid->niib2);
	  jj = is2; /* Choose outflow as default */
	}
	for (n=0; n<NVAR; n++) {
	  bval->uiib2[i][j][n] = agrid->u[i][jj][n];
	}
#ifdef MHD
	bval->bxiib2[i][j] = agrid->bx[i][jj];
	bval->byiib2[i][j] = agrid->by[i][jj];
#endif /* MHD */
	bval->uiib2[i][j][2] = agrid->u[i][jj][2]*reflect_s;
      }
    }
  }

/*----------------------- Continue with Outer 2 boundary ---------------------*/

  if (agrid->noib2 == 1) {reflect_s = -1.0;} else {reflect_s = 1.0;}

/* Do noib2=3 (inflow) */

  if (agrid->noib2 == 3) {
    for (i=is1-4; i<=ie1+4; i++){
      for (j=0; j<4; j++) {
	for (n=0; n<NVAR; n++) {
	  bval->uoib2[i][j][n] = agrid->boundary_values.uoib2[i][j][n];
	}
#ifdef MHD
	bval->bxoib2[i][j] = agrid->boundary_values.bxoib2[i][j];
	bval->byoib2[i][j] = agrid->boundary_values.byoib2[i][j];
#endif /* MHD */
      }
    }

/* Do noib2=1 (reflection), 2 (flow-out), 4 (periodic) by adjusting i-index   */
/* Note Bx ghost zones start at ie+2 rather than ie+1 */

  } else {
    for(i=is1-4; i<=ie1+4; i++){
      for (j=0; j<4; j++) {
	switch(agrid->noib2){
	case 1:
	  jj = ie2-j;  break;
	case 2:
	  jj = ie2;  break;
	case 4:
	  jj = is2+j;  break;
	default:
	  fprintf(stderr,"[set_bval_arrays]: agrid->noib2 = %d unknown\n",
		  agrid->noib2);
	  jj = ie2; /* Choose outflow as default */
	}
	for (n=0; n<NVAR; n++) {
	  bval->uoib2[i][j][n] = agrid->u[i][jj][n];
	}
#ifdef MHD
	bval->bxoib2[i][j] = agrid->bx[i][jj];
	bval->byoib2[i][j] = agrid->by[i][jj+1];
#endif /* MHD */
	bval->uoib2[i][j][2] = agrid->u[i][jj][2]*reflect_s;
      }
    }
  }

}


#else /* THREE_D */
#error : [set_bval_arrays] not yet implemented for 3D
#endif

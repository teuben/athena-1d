/*==============================================================================
 *                             Function DATA_OUTPUT 
 * Outputs data in input grid_block structure in various formats
 * On input, *_flag = 1 forces dump of type *
 *============================================================================*/
#include <stdio.h>
#include "athena.def"
#include "athena.h"
void data_output(struct grid_block *agrid, int hst_flag, int hdf_flag,
  int bin_flag)
{
#include "prototypes.h"

/* History Dump */

   if (agrid->time >= (agrid->t_hst + agrid->dt_hst)) {
      agrid->t_hst = agrid->t_hst + agrid->dt_hst;
      hst_flag = 1;
   }
   if (hst_flag == 1) printd(agrid);

/* Unformatted (binary) Dump */

   if (agrid->time >= (agrid->t_bin + agrid->dt_bin)) {
      agrid->t_bin = agrid->t_bin + agrid->dt_bin;
      bin_flag = 1;
   }
   if (bin_flag == 1) binary_dump(agrid);

}

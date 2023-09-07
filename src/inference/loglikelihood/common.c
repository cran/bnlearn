#include "../../include/rcore.h"
#include "../../fitted/fitted.h"
#include "../../core/data.table.h"

/* check whether the data are complete for all the local distributions of
 * interest. */
bool check_locally_incomplete_data(fitted_bn bn, meta m, bool debugging) {

bool early_return = FALSE;

  for (int i = 0; i < m.ncols; i++) {

    if (!m.flag[i].fixed)
      continue;

    /* for the data to be locally complete, the node itself... */
    if (!m.flag[i].complete) {

      early_return  = TRUE;
      goto incomplete_data;

    }/*THEN*/

    /* ... and its parents must be complete. */
    for (int j = 0; j < bn.ldists[i].nparents; j++)
       if (!m.flag[bn.ldists[i].parents[j]].complete) {

         early_return = TRUE;
         goto incomplete_data;

       }/*THEN*/

incomplete_data:
    if (early_return) {

      if (debugging)
        Rprintf("* incomplete data for node %s, the log-likelihood is NA.\n",
          bn.labels[i]);

      return TRUE;

    }/*THEN*/

  }/*FOR*/

  return FALSE;

}/*CHECK_LOCALLY_INCOMPLETE_DATA*/


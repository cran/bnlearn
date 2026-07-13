#include "../include/rcore.h"
#include "../minimal/common.h"
#include "allocations.h"
#include "data.table.h"
#include "moments.h"

/* -------------------- metadata ----------------------------------------- */

void meta_init_flags(meta *m, int offset, SEXP complete, SEXP fixed) {

  if (complete != R_NilValue) {

    int *cc = LOGICAL(complete);

    for (int i = 0; i < (*m).ncols; i++)
      (*m).flag[i + offset].complete = cc[i];

  }/*THEN*/

  if (fixed != R_NilValue) {

    int *ff = INTEGER(fixed);

    for (int i = 0; i < length(fixed); i++)
      (*m).flag[i + offset].fixed = (ff[i] > 0);

  }/*THEN*/

}/*META_INIT_FLAGS*/

void meta_copy_names(meta *m, int offset, SEXP df) {

SEXP names = getAttrib(df, R_NamesSymbol);

  if (!(*m).names)
    (*m).names = Calloc1D((*m).ncols, sizeof(char *));

  for (int i = 0; i + offset < (*m).ncols; i++)
    (*m).names[i + offset] = CHAR(STRING_ELT(names, i));

}/*META_COPY_NAMES*/

void meta_drop_flagged(meta *src, meta *dest) {

int i = 0, j = 0;

  if (!(*dest).names && (*src).names)
    (*dest).names = Calloc1D((*src).ncols, sizeof(char *));

  for (i = 0; i < (*src).ncols; i++) {

    if ((*src).flag[i].drop)
      continue;

    if ((*src).names)
      (*dest).names[j] = (*src).names[i];
    if ((*src).flag)
      (*dest).flag[j] = (*src).flag[i];

    /* the source and the destination are different structs, the columns belong
     * with the former and not with the latter. */
    if (src != dest)
      (*dest).flag[j].own = FALSE;

    j++;

  }/*FOR*/

  (*dest).nobs = (*src).nobs;
  (*dest).ncols = j;

}/*META_DROP_FLAGGED*/

void meta_subset_columns(meta *src, meta *dest, int *ids, int nids) {

int i = 0;

  if (!(*dest).names && (*src).names)
    (*dest).names = Calloc1D((*src).ncols, sizeof(char *));

  for (i = 0; i < nids; i++) {

    if ((*src).names)
      (*dest).names[i] = (*src).names[ids[i]];
    if ((*src).flag)
      (*dest).flag[i] = (*src).flag[ids[i]];

    /* the source and the destination are different structs, the columns belong
     * with the former and not with the latter. */
    if (src != dest)
      (*dest).flag[i].own = FALSE;

  }/*FOR*/

  (*dest).nobs = (*src).nobs;
  (*dest).ncols = nids;

}/*META_SUBSET_COLUMNS*/

void meta_copy(meta *src, meta *dest) {

int i = 0;
bool own = FALSE;

  /* copy all the flags that describe the contents and the types of the columns,
   * but preserve ownership of the memory in the destination struct. */
  for (i = 0; i < (*src).ncols; i++) {

    own = (*dest).flag[i].own;
    (*dest).flag[i] = (*src).flag[i];
    (*dest).flag[i].own = own;

  }/*FOR*/

  (*dest).nobs = (*src).nobs;
  (*dest).ncols = (*src).ncols;

}/*META_COPY*/

void print_meta(meta *m, int i) {

  Rprintf("%10s", (*m).names ? (*m).names[i] : "");
  Rprintf(" [%s%s%s%s%s%s]",
    ((*m).flag[i].own ? "O" : "P"),
    ((*m).flag[i].discrete ? "D" : " "),
    ((*m).flag[i].continuous ? "C" : " "),
    ((*m).flag[i].complete ? "C" : " "),
    ((*m).flag[i].fixed ? "F" : " "),
    ((*m).flag[i].drop ? "D" : " "));

}/*PRINT_META*/

void FreeMETA(meta *m) {

  Free1D((*m).flag);
  Free1D((*m).names);

}/*FREEMETA*/

/* -------------------- tabular data structure --------------------------- */

/* create a unified data table from an R data frame. */
tabular tabular_from_SEXP(SEXP df, int doffset, int coffset) {

int i = 0, j = 0, k = 0, nc = 0, nd = 0, nrows = 0, ncols = length(df);
SEXP *vec = NULL;
tabular dt = { 0 };

  vec = Calloc1D(ncols, sizeof(SEXP));

  /* count how many columns of each type are present in the data frame. */
  for (i = 0; i < ncols; i++) {

    vec[i] = VECTOR_ELT(df, i);

    if (TYPEOF(vec[i]) == REALSXP)
      nc++;
    else
      nd++;

  }/*FOR*/

  /* initialize, taking care of the corner case of zero-columns data frames. */
  if (ncols > 0)
    nrows = length(VECTOR_ELT(df, 0));
  else {

    SEXP dim = getAttrib(df, R_RowNamesSymbol);
    if (isInteger(dim) && length(dim) == 2 && INTEGER(dim)[0] == NA_INTEGER)
      nrows = abs(INTEGER(dim)[1]);
    else
      nrows = length(dim);

  }/*ELSE*/

  dt = empty_tabular(nrows, nd + doffset, nc + coffset);

  /* reserve the leading continuous and discrete slots (the offsets). */
  for (i = 0; i < coffset; i++) {

    dt.map[i] = i;
    dt.m.flag[i].continuous = TRUE;

  }/*FOR*/

  for (i = 0; i < doffset; i++) {

    dt.map[i + coffset] = i;
    dt.m.flag[i + coffset].discrete = TRUE;

  }/*FOR*/

  for (i = 0, j = coffset, k = doffset; i < ncols; i++) {

    switch(TYPEOF(vec[i])) {

      case REALSXP:
        dt.ccol[j] = REAL(vec[i]);
        dt.map[i + coffset + doffset] = j;
        dt.m.flag[i + coffset + doffset].continuous = TRUE;
        j++;
        break;

      case INTSXP:
        dt.dcol[k] = INTEGER(vec[i]);
        /* save the numbers of levels as well for factors, to iterate over
         * contingency tables, compute configurations and degrees of freedom. */
        dt.nlvl[k] = NLEVELS(vec[i]);
        dt.map[i + coffset + doffset] = k;
        dt.m.flag[i + coffset + doffset].discrete = TRUE;
        k++;
        break;

      default:
        error("this SEXP type is not handled in tabular_from_SEXP().");

    }/*SWITCH*/

  }/*FOR*/

  Free1D(vec);

  return dt;

}/*TABULAR_FROM_SEXP*/

/* barebones tabular object (no data or allocated memory, minimal metadata). */
tabular empty_tabular(int nobs, int dcols, int ccols) {

int i = 0;
tabular dt = { 0 };

  dt.m.nobs = nobs;
  dt.m.ncols = dcols + ccols;
  dt.nccols = ccols;
  dt.ccol = Calloc1D(ccols, sizeof(double *));
  dt.ndcols = dcols;
  dt.dcol = Calloc1D(dcols, sizeof(int *));
  dt.nlvl = Calloc1D(dcols, sizeof(int));
  dt.map = Calloc1D(dcols + ccols, sizeof(int));

  dt.m.flag = Calloc1D(dcols + ccols, sizeof(flags));
  for (i = 0; i < dcols + ccols; i++)
    dt.m.flag[i] = (flags){ .complete = TRUE, .own = FALSE };

  return dt;

}/*EMPTY_TABULAR*/

/* empty tabular object (with memory allocated to hold data). */
tabular new_tabular(int nobs, int dcols, int ccols) {

int j = 0;
tabular dt = empty_tabular(nobs, dcols, ccols);

  for (j = 0; j < dcols; j++)
    dt.dcol[j] = Calloc1D(nobs, sizeof(int));
  for (j = 0; j < ccols; j++)
    dt.ccol[j] = Calloc1D(nobs, sizeof(double));

  /* the columns belong with the struct, flag them so that they are freed. */
  for (j = 0; j < dcols + ccols; j++)
    dt.m.flag[j].own = TRUE;

  /* map the columns to the pointers, sequentially, starting from the discrete
   * variables. */
  for (j = 0; j < dcols; j++)
    dt.map[j] = j;
  for (j = 0; j < ccols; j++)
    dt.map[dcols + j] = j;

  /* set the flags for the data types to match the map. */
  for (j = 0; j < dcols; j++) {

    dt.m.flag[j].discrete = TRUE;
    dt.m.flag[j].continuous = FALSE;

  }/*FOR*/
  for (j = 0; j < ccols; j++) {

    dt.m.flag[dcols + j].discrete = FALSE;
    dt.m.flag[dcols + j].continuous = TRUE;

  }/*FOR*/

  return dt;

}/*NEW_TABULAR*/

/* cache the means of the continuous columns. */
void tabular_cache_means(tabular *dt, int offset) {

  (*dt).mean = Calloc1D((*dt).nccols, sizeof(double));
  c_meanvec((*dt).ccol, (*dt).mean, (*dt).m.nobs, (*dt).nccols, offset);

}/*TABULAR_CACHE_MEANS*/

/* print tabular objects in C code. */
void print_tabular(tabular dt) {

int i = 0;

  Rprintf("tabular: %dx%d (%d discrete, %d continuous) \n",
    dt.m.nobs, dt.m.ncols, dt.ndcols, dt.nccols);

  for (i = 0; i < dt.m.ncols; i++) {

    print_meta(&(dt.m), i);

    if (dt.m.flag[i].discrete)
      Rprintf(" levels: %d", dt.nlvl[dt.map[i]]);
    else if (dt.m.flag[i].continuous) {

      Rprintf("@%p", (void *)dt.ccol[dt.map[i]]);
      if (dt.mean)
        Rprintf(" mean: %lf", dt.mean[dt.map[i]]);

    }/*THEN*/

    Rprintf("\n");

  }/*FOR*/

}/*PRINT_TABULAR*/

/* remove flagged columns from a tabular object. */
void tabular_drop_flagged(tabular *dt, tabular *copy) {

int i = 0, j = 0, k = 0;

  /* subset the data. */
  for (i = 0; i < (*dt).m.ncols; i++) {

    /* do not copy columns marked as to be dropped. */
    if ((*dt).m.flag[i].drop)
      continue;

    if ((*dt).m.flag[i].discrete) {

      (*copy).dcol[j] = (*dt).dcol[(*dt).map[i]];
      (*copy).nlvl[j] = (*dt).nlvl[(*dt).map[i]];
      (*copy).map[j + k] = j;
      j++;

    }/*THEN*/
    else if ((*dt).m.flag[i].continuous) {

      (*copy).ccol[k] = (*dt).ccol[(*dt).map[i]];
      if ((*dt).mean && (*copy).mean)
        (*copy).mean[k] = (*dt).mean[(*dt).map[i]];
      (*copy).map[j + k] = k;
      k++;

    }/*THEN*/

  }/*FOR*/

  /* update the columns counts. */
  (*copy).ndcols = j;
  (*copy).nccols = k;

  /* subset the metadata. */
  meta_drop_flagged(&((*dt).m), &((*copy).m));

}/*TABULAR_DROP_FLAGGED*/

/* copy a subset of columns to a second tabular object. */
void tabular_subset_columns(tabular *dt, tabular *copy, int *ids, int nids) {

int i = 0, j = 0, k = 0;

  /* subset the data. */
  for (i = 0; i < nids; i++) {

    if ((*dt).m.flag[ids[i]].discrete) {

      (*copy).dcol[j] = (*dt).dcol[(*dt).map[ids[i]]];
      (*copy).nlvl[j] = (*dt).nlvl[(*dt).map[ids[i]]];
      (*copy).map[j + k] = j;
      j++;

    }/*THEN*/
    else if ((*dt).m.flag[ids[i]].continuous) {

      (*copy).ccol[k] = (*dt).ccol[(*dt).map[ids[i]]];
      if ((*dt).mean && (*copy).mean)
        (*copy).mean[k] = (*dt).mean[(*dt).map[ids[i]]];
      (*copy).map[j + k] = k;
      k++;

    }/*THEN*/

  }/*FOR*/

  /* update the columns counts. */
  (*copy).ndcols = j;
  (*copy).nccols = k;

  /* subset the metadata. */
  meta_subset_columns(&((*dt).m), &((*copy).m), ids, nids);

}/*TABULAR_SUBSET_COLUMNS*/

/* flag incomplete cases looking at a subset of columns. */
void tabular_incomplete_cases_range(tabular *dt, bool *indicator, int dcol_start,
    int dcol_end, int ccol_start, int ccol_end) {

int i = 0, j = 0;

  for (j = dcol_start; j <= dcol_end; j++)
    for (i = 0; i < (*dt).m.nobs; i++)
      if ((*dt).dcol[j][i] == NA_INTEGER)
        indicator[i] = TRUE;

  for (j = ccol_start; j <= ccol_end; j++)
    for (i = 0; i < (*dt).m.nobs; i++)
      if (ISNAN((*dt).ccol[j][i]))
        indicator[i] = TRUE;

}/*TABULAR_INCOMPLETE_CASES_RANGE*/

/* flag incomplete cases looking at all columns. */
void tabular_incomplete_cases(tabular *dt, bool *indicator, int doffset,
    int coffset) {

  tabular_incomplete_cases_range(dt, indicator, doffset, (*dt).ndcols - 1,
    coffset, (*dt).nccols - 1);

}/*TABULAR_INCOMPLETE_CASES*/

/* copy a subset of observations to a second tabular object. */
void tabular_subsample_by_logical(tabular *dt, tabular *copy, bool *indicators,
    int doffset, int coffset) {

int j = 0, k = 0, l = 0;

  if ((*dt).m.ncols == 0) {

    for (l = 0, k = 0; l < (*dt).m.nobs; l++)
      if (!indicators[l])
        k++;

  }/*THEN*/
  else {

    for (j = coffset; j < (*dt).nccols; j++)
      for (l = 0, k = 0; l < (*dt).m.nobs; l++)
        if (!indicators[l])
          (*copy).ccol[j][k++] = (*dt).ccol[j][l];

    for (j = doffset; j < (*dt).ndcols; j++)
      for (l = 0, k = 0; l < (*dt).m.nobs; l++)
        if (!indicators[l])
          (*copy).dcol[j][k++] = (*dt).dcol[j][l];

  }/*ELSE*/

  meta_copy(&((*dt).m), &((*copy).m));
  (*copy).m.nobs = k;
  (*copy).ndcols = (*dt).ndcols;
  (*copy).nccols = (*dt).nccols;

  for (j = 0; j < (*dt).ndcols; j++)
    (*copy).nlvl[j] = (*dt).nlvl[j];

  for (j = 0; j < (*dt).m.ncols; j++)
    (*copy).map[j] = (*dt).map[j];

  if ((*dt).mean && (*copy).mean)
    for (j = 0; j < (*dt).nccols; j++)
      (*copy).mean[j] = (*dt).mean[j];

  if ((*dt).m.names && (*copy).m.names)
    for (j = 0; j < (*dt).m.ncols; j++)
      (*copy).m.names[j] = (*dt).m.names[j];

}/*TABULAR_SUBSAMPLE_BY_LOGICAL*/

/* free all memory used by a tabular object. */
void FreeTAB(tabular dt) {

int j = 0;

  /* free the columns that belong with the struct, and leave the rest alone. */
  for (j = 0; j < dt.m.ncols; j++) {

    if (!dt.m.flag[j].own)
      continue;

    if (dt.m.flag[j].discrete)
      Free1D(dt.dcol[dt.map[j]]);
    else if (dt.m.flag[j].continuous)
      Free1D(dt.ccol[dt.map[j]]);

  }/*FOR*/

  /* free the column pointers, unconditionally. */
  Free1D(dt.ccol);
  Free1D(dt.dcol);

  /* free the cached means, the numbers of levels, the column map and, finally,
   * the meta data. */
  Free1D(dt.mean);
  Free1D(dt.nlvl);
  Free1D(dt.map);
  FreeMETA(&(dt.m));

}/*FREETAB*/


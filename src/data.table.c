#include "include/rcore.h"
#include "include/data.table.h"
#include "include/covariance.h"

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
  Rprintf(" [%s%s%s%s%s]",
    ((*m).flag[i].own ? "O" : "P"),
    ((*m).flag[i].discrete ? "D" : " "),
    ((*m).flag[i].gaussian ? "G" : " "),
    ((*m).flag[i].complete ? "C" : " "),
    ((*m).flag[i].fixed ? "F" : " "),
    ((*m).flag[i].drop ? "D" : " "));

}/*PRINT_META*/

void FreeMETA(meta *m) {

  Free1D((*m).flag);
  Free1D((*m).names);

}/*FREEMETA*/

/* -------------------- discrete data ------------------------------------ */

ddata ddata_from_SEXP(SEXP df, int offset) {

int i = 0, ncols = length(df);
SEXP col;
ddata dt = { 0 };

  /* initialize the object, taking care of the corner case of zero-columns data
   * frames. */
  if (ncols == 0)
    dt = empty_ddata(0, ncols + offset);
  else
    dt = empty_ddata(length(VECTOR_ELT(df, 0)), ncols + offset);

  for (i = 0; i < ncols; i++) {

    col = VECTOR_ELT(df, i);
    dt.col[i + offset] = INTEGER(col);
    dt.nlvl[i + offset] = NLEVELS(col);

  }/*FOR*/

  return dt;

}/*DDATA_FROM_SEXP*/

ddata empty_ddata(int nobs, int ncols) {

int i = 0;
ddata dt = { 0 };

  dt.m.nobs = nobs;
  dt.m.ncols = ncols;
  dt.col = Calloc1D(ncols, sizeof(int *));
  dt.nlvl = Calloc1D(ncols, sizeof(int));

  dt.m.flag = Calloc1D(ncols, sizeof(flags));
  for (i = 0; i < ncols; i++)
    dt.m.flag[i] = (flags){ .discrete = TRUE, .complete = TRUE, .own = FALSE };

  return dt;

}/*EMPTY_DDATA*/

void print_ddata(ddata dt) {

int i = 0;

  Rprintf("ddata: %dx%d\n", dt.m.nobs, dt.m.ncols);

  for (i = 0; i < dt.m.ncols; i++) {

    print_meta(&(dt.m), i);
    Rprintf("@%p", (void *)dt.col[i]);
      Rprintf(" levels: %d", dt.nlvl[i]);

    Rprintf("\n");

  }/*FOR*/

}/*PRINT_DDATA*/

void ddata_drop_flagged(ddata *dt, ddata *copy) {

int i = 0, j = 0;

  /* subset the data. */
  for (i = 0; i < (*dt).m.ncols; i++) {

    /* do not copy columns marked as to be dropped. */
    if ((*dt).m.flag[i].drop)
      continue;

    (*copy).col[j] = (*dt).col[i];
    (*copy).nlvl[j] = (*dt).nlvl[i];
    j++;

  }/*FOR*/

  /* subset the metadata. */
  meta_drop_flagged(&((*dt).m), &((*copy).m));

}/*DDATA_DROP_FLAGGED*/

void ddata_subset_columns(ddata *dt, ddata *copy, int *ids, int nids) {

int i = 0;

  /* subset the data. */
  for (i = 0; i < nids; i++) {

    (*copy).col[i] = (*dt).col[ids[i]];
    (*copy).nlvl[i] = (*dt).nlvl[ids[i]];

  }/*FOR*/

  /* subset the metadata. */
  meta_subset_columns(&((*dt).m), &((*copy).m), ids, nids);

}/*DDATA_SUBSET_COLUMNS*/

/* free a discrete data table. */
void FreeDDT(ddata dt) {

int i = 0;

  /* free the columns that belong with the struct, and leave the rest alone. */
  for (i = 0; i < dt.m.ncols; i++)
    if (dt.m.flag[i].own)
      Free1D(dt.col[i]);
  /* free the column pointers, unconditionally. */
  Free1D(dt.col);

  /* free the numbers of levels and, finally, the meta data. */
  Free1D(dt.nlvl);
  FreeMETA(&(dt.m));

}/*FREEDDT*/

/* -------------------- Gaussian data ------------------------------------ */

/* create a data table from Gaussian data. */
gdata gdata_from_SEXP(SEXP df, int offset) {

int i = 0, ncols = length(df);
gdata dt = { 0 };

  /* initialize the object, taking care of the corner case of zero-columns data
   * frames. */
  if (ncols == 0)
    dt = empty_gdata(0, ncols + offset);
  else
    dt = empty_gdata(length(VECTOR_ELT(df, 0)), ncols + offset);

  for (i = 0; i < ncols; i++)
    dt.col[i + offset] = REAL(VECTOR_ELT(df, i));

  return dt;

}/*GDATA_FROM_SEXP*/

gdata new_gdata(int nobs, int ncols) {

gdata dt = empty_gdata(nobs, ncols);

  for (int j = 0; j < ncols; j++) {

    dt.col[j] = Calloc1D(nobs, sizeof(double));
    /* the columns belong with the struct, flag them so that they are freed. */
    dt.m.flag[j].own = TRUE;

  }/*FOR*/

  return dt;

}/*NEW_GDATA*/

gdata empty_gdata(int nobs, int ncols) {

int i = 0;
gdata dt = { 0 };

  dt.m.nobs = nobs;
  dt.m.ncols = ncols;
  dt.col = Calloc1D(ncols, sizeof(double *));

  dt.m.flag = Calloc1D(ncols, sizeof(flags));
  for (i = 0; i < ncols; i++)
    dt.m.flag[i] = (flags){ .gaussian = TRUE, .complete = TRUE, .own = FALSE };

  return dt;

}/*EMPTY_GDATA*/

/* cache the means of the continuous columns. */
void gdata_cache_means(gdata *dt, int offset) {

  (*dt).mean = Calloc1D((*dt).m.ncols, sizeof(double));
  c_meanvec((*dt).col, (*dt).mean, (*dt).m.nobs, (*dt).m.ncols, offset);

}/*CACHE_MEANS*/

void print_gdata(gdata dt) {

int i = 0;

  Rprintf("gdata: %dx%d\n", dt.m.nobs, dt.m.ncols);

  for (i = 0; i < dt.m.ncols; i++) {

    print_meta(&(dt.m), i);
    Rprintf("@%p", (void *)dt.col[i]);
    if (dt.mean)
      Rprintf(" mean: %lf", dt.mean[i]);

    Rprintf("\n");

  }/*FOR*/

}/*PRINT_GDATA*/

void gdata_drop_flagged(gdata *dt, gdata *copy) {

int i = 0, k = 0;

  /* subset the data. */
  for (i = 0; i < (*dt).m.ncols; i++) {

    /* do not copy columns marked as to be dropped. */
    if ((*dt).m.flag[i].drop)
      continue;

    (*copy).col[k] = (*dt).col[i];
    if ((*dt).mean && (*copy).mean)
      (*copy).mean[k] = (*dt).mean[i];
    k++;

  }/*FOR*/

  /* subset the metadata. */
  meta_drop_flagged(&((*dt).m), &((*copy).m));

}/*GDATA_DROP_FLAGGED*/

void gdata_subset_columns(gdata *dt, gdata *copy, int *ids, int nids) {

int i = 0;

  /* subset the data. */
  for (i = 0; i < nids; i++) {

    (*copy).col[i] = (*dt).col[ids[i]];
    if ((*dt).mean && (*copy).mean)
      (*copy).mean[i] = (*dt).mean[ids[i]];

  }/*FOR*/

  /* subset the metadata. */
  meta_subset_columns(&((*dt).m), &((*copy).m), ids, nids);

}/*GDATA_SUBSET_COLUMNS*/

void gdata_incomplete_cases_range(gdata *dt, bool *indicator, int col_start,
    int col_end) {

int i = 0, j = 0;

  for (i = 0; i < (*dt).m.nobs; i++)
    for (j = col_start; j <= col_end; j++)
      if (ISNAN((*dt).col[j][i])) {

        indicator[i] = TRUE;
        break;

      }/*THEN*/

}/*GDATA_INCOMPLETE_CASES*/

void gdata_subsample_by_logical(gdata *dt, gdata *copy, bool *indicators,
    int offset) {

int i = 0, j = 0, k = 0;

  for (j = offset; j < (*dt).m.ncols; j++)
    for (i = 0, k = 0; i < (*dt).m.nobs; i++)
      if (!indicators[i])
        (*copy).col[j][k++] = (*dt).col[j][i];

  meta_copy(&((*dt).m), &((*copy).m));
  (*copy).m.nobs = k;

  if ((*dt).m.names && (*copy).m.names)
    for (j = 0; j < (*dt).m.ncols; j++)
      (*copy).m.names[j] = (*dt).m.names[j];

}/*GDATA_SUBSAMPLE_BY_LOGICAL*/

/* free a Gaussian data table. */
void FreeGDT(gdata dt) {

int i = 0;

  /* free the columns that belong with the struct, and leave the rest alone. */
  for (i = 0; i < dt.m.ncols; i++)
    if (dt.m.flag[i].own)
      Free1D(dt.col[i]);
  /* free the column pointers, unconditionally. */
  Free1D(dt.col);

  /* free the cached means and, finally, the meta data. */
  Free1D(dt.mean);
  FreeMETA(&(dt.m));

}/*FREEGDT*/

/* -------------------- conditional Gaussian data ------------------------ */

/* create a data table from conditional Gaussian data. */
cgdata cgdata_from_SEXP(SEXP df, int doffset, int goffset) {

int i = 0, j = 0, k = 0, nc = 0, nd = 0, ncols = length(df);
SEXP *vec = NULL;
cgdata dt = { 0 };

  vec = Calloc1D(ncols, sizeof(SEXP));

  /* count how many columns of each type are present in the data frame. */
  for (i = 0; i < ncols; i++) {

    vec[i] = VECTOR_ELT(df, i);

    if (TYPEOF(vec[i]) == REALSXP)
      nc++;
    else
      nd++;

  }/*FOR*/

  /* initialize the object, taking care of the corner case of zero-columns data
   * frames. */
  if (ncols == 0)
    dt = empty_cgdata(0, nd + doffset, nc + goffset);
  else
    dt = empty_cgdata(length(vec[0]), nd + doffset, nc + goffset);

  for (i = 0; i < goffset; i++) {

    dt.map[i] = i;
    dt.m.flag[i].gaussian = TRUE;

  }/*FOR*/

  for (i = 0; i < doffset; i++) {

    dt.map[i + goffset] = i;
    dt.m.flag[i + goffset].discrete = TRUE;

  }/*FOR*/

  for (i = 0, j = goffset, k = doffset; i < ncols; i++) {

    switch(TYPEOF(vec[i])) {

      case REALSXP:
        dt.gcol[j] = REAL(vec[i]);
        dt.map[i + goffset + doffset] = j;
        dt.m.flag[i + goffset + doffset].gaussian = TRUE;
        j++;
        break;

      case INTSXP:
        dt.dcol[k] = INTEGER(vec[i]);
        /* save the numbers of levels as well for factors, to iterate over
         * contingency tables, compute configurations and degrees of freedom. */
        dt.nlvl[k] = NLEVELS(vec[i]);
        dt.map[i + goffset + doffset] = k;
        dt.m.flag[i + goffset + doffset].discrete = TRUE;
        k++;
        break;

      default:
        error("this SEXP type is not handled in data_table_from_SEXP().");

    }/*SWITCH*/

  }/*FOR*/

  Free1D(vec);

  return dt;

}/*CGDATA_FROM_SEXP*/

cgdata new_cgdata(int nobs, int dcols, int gcols) {

int j = 0;
cgdata dt = empty_cgdata(nobs, dcols, gcols);

  for (j = 0; j < dcols; j++)
    dt.dcol[j] = Calloc1D(nobs, sizeof(int));
  for (j = 0; j < gcols; j++)
    dt.gcol[j] = Calloc1D(nobs, sizeof(double));

  /* the columns belong with the struct, flag them so that they are freed. */
  for (j = 0; j < dcols + gcols; j++)
    dt.m.flag[j].own = TRUE;

  /* map the columns to the pointes, sequentially, starting from the discrete
   * variables. */
  for (j = 0; j < dcols; j++)
    dt.map[j] = j;
  for (j = 0; j < gcols; j++)
    dt.map[dcols + j] = j;

  /* set the flags for the data types to match the map. */
  for (j = 0; j < dcols; j++) {

    dt.m.flag[j].discrete = TRUE;
    dt.m.flag[j].gaussian = FALSE;

  }/*FOR*/
  for (j = 0; j < gcols; j++) {

    dt.m.flag[dcols + j].discrete = FALSE;
    dt.m.flag[dcols + j].gaussian = TRUE;

  }/*FOR*/

  return dt;

}/*NEW_CGDATA*/

cgdata empty_cgdata(int nobs, int dcols, int gcols) {

int i = 0;
cgdata dt = { 0 };

  dt.m.nobs = nobs;
  dt.m.ncols = dcols + gcols;
  dt.ngcols = gcols;
  dt.gcol = Calloc1D(gcols, sizeof(double *));
  dt.ndcols = dcols;
  dt.dcol = Calloc1D(dcols, sizeof(int *));
  dt.nlvl = Calloc1D(dcols, sizeof(int));
  dt.map = Calloc1D(dcols + gcols, sizeof(int));

  dt.m.flag = Calloc1D(dcols + gcols, sizeof(flags));
  for (i = 0; i < dcols + gcols; i++)
    dt.m.flag[i] = (flags){ .complete = TRUE, .own = FALSE };

  return dt;

}/*EMPTY_CGDATA*/

void print_cgdata(cgdata dt) {

int i = 0;

  Rprintf("cgdata: %dx%d (%d discrete, %d continuous) \n",
    dt.m.nobs, dt.m.ncols, dt.ndcols, dt.ngcols);

  for (i = 0; i < dt.m.ncols; i++) {

    print_meta(&(dt.m), i);

    if (dt.m.flag[i].discrete)
      Rprintf(" levels: %d", dt.nlvl[dt.map[i]]);
    else
      Rprintf("@%p", (void *)dt.gcol[dt.map[i]]);

    Rprintf("\n");

  }/*FOR*/

}/*PRINT_CGDATA*/

void cgdata_drop_flagged(cgdata *dt, cgdata *copy) {

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
    else if ((*dt).m.flag[i].gaussian) {

      (*copy).gcol[k] = (*dt).gcol[(*dt).map[i]];
      (*copy).map[j + k] = k;
      k++;

    }/*THEN*/

  }/*FOR*/

  /* update the columns counts. */
  (*copy).ndcols = j;
  (*copy).ngcols = k;

  /* subset the metadata. */
  meta_drop_flagged(&((*dt).m), &((*copy).m));

}/*CGDATA_DROP_FLAGGED*/

void cgdata_subset_columns(cgdata *dt, cgdata *copy, int *ids, int nids) {

int i = 0, j = 0, k = 0;

  /* subset the data. */
  for (i = 0; i < nids; i++) {

    if ((*dt).m.flag[ids[i]].discrete) {

      (*copy).dcol[j] = (*dt).dcol[(*dt).map[ids[i]]];
      (*copy).nlvl[j] = (*dt).nlvl[(*dt).map[ids[i]]];
      (*copy).map[j + k] = j;
      j++;

    }/*THEN*/
    else if ((*dt).m.flag[ids[i]].gaussian) {

      (*copy).gcol[k] = (*dt).gcol[(*dt).map[ids[i]]];
      (*copy).map[j + k] = k;
      k++;

    }/*THEN*/

  }/*FOR*/

  /* update the columns counts. */
  (*copy).ndcols = j;
  (*copy).ngcols = k;

  /* subset the metadata. */
  meta_subset_columns(&((*dt).m), &((*copy).m), ids, nids);

}/*CGDATA_SUBSET_COLUMNS*/

void cgdata_incomplete_cases(cgdata *dt, bool *indicator, int doffset, int goffset) {

int i = 0, j = 0;

  for (i = 0; i < (*dt).m.nobs; i++) {

    for (j = doffset; j < (*dt).ndcols; j++)
      if ((*dt).dcol[j][i] == NA_INTEGER) {

        indicator[i] = TRUE;
        continue;

      }/*THEN*/

    for (j = goffset; j < (*dt).ngcols; j++)
      if (ISNAN((*dt).gcol[j][i])) {

        indicator[i] = TRUE;
        continue;

      }/*THEN*/

  }/*FOR*/

}/*CGDATA_INCOMPLETE_CASES*/

void cgdata_subsample_by_logical(cgdata *dt, cgdata *copy, bool *indicators,
    int doffset, int goffset) {

int j = 0, k = 0, l = 0;

  for (j = doffset; j < (*dt).ngcols; j++)
    for (l = 0, k = 0; l < (*dt).m.nobs; l++)
      if (!indicators[l])
        (*copy).gcol[j][k++] = (*dt).gcol[j][l];

  for (j = goffset; j < (*dt).ndcols; j++)
    for (l = 0, k = 0; l < (*dt).m.nobs; l++)
      if (!indicators[l])
        (*copy).dcol[j][k++] = (*dt).dcol[j][l];

  meta_copy(&((*dt).m), &((*copy).m));
  (*copy).m.nobs = k;
  (*copy).ndcols = (*dt).ndcols;
  (*copy).ngcols = (*dt).ngcols;

  for (j = 0; j < (*dt).ndcols; j++)
    (*copy).nlvl[j] = (*dt).nlvl[j];

  for (j = 0; j < (*dt).m.ncols; j++)
    (*copy).map[j] = (*dt).map[j];

  if ((*dt).m.names && (*copy).m.names)
    for (j = 0; j < (*dt).m.ncols; j++)
      (*copy).m.names[j] = (*dt).m.names[j];

}/*CGDATA_SUBSAMPLE_BY_LOGICAL*/

void FreeCGDT(cgdata dt) {

int j = 0;

  /* free the columns that belong with the struct, and leave the rest alone. */
  for (j = 0; j < dt.m.ncols; j++) {

    if (!dt.m.flag[j].own)
      continue;

    if (dt.m.flag[j].discrete)
      Free1D(dt.dcol[dt.map[j]]);
    else if (dt.m.flag[j].gaussian)
      Free1D(dt.gcol[dt.map[j]]);

  }/*FOR*/

  /* free the column pointers, unconditionally. */
  Free1D(dt.gcol);
  Free1D(dt.dcol);

  /* free the numbers of levels, the column map and, finally, the meta data. */
  Free1D(dt.nlvl);
  Free1D(dt.map);
  FreeMETA(&(dt.m));

}/*FREECGDT*/


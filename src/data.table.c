#include "include/rcore.h"
#include "include/data.table.h"
#include "include/covariance.h"

/* -------------------- metadata ----------------------------------------- */

void meta_init_flags(meta *m, int offset, SEXP complete, SEXP fixed) {

  /* allocate the flags if not allocated already (as is the case for cgdata). */
  if (!(*m).flag)
    (*m).flag = Calloc1D((*m).ncols, sizeof(flags));

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
  if (!(*dest).flag && (*src).flag)
    (*dest).flag = Calloc1D((*src).ncols, sizeof(flags));

  for (i = 0; i < (*src).ncols; i++) {

    if ((*src).flag[i].drop)
      continue;

    if ((*src).names)
      (*dest).names[j] = (*src).names[i];
    if ((*src).flag)
      (*dest).flag[j] = (*src).flag[i];
    j++;

  }/*FOR*/

  (*dest).nobs = (*src).nobs;
  (*dest).ncols = j;

}/*META_DROP_FLAGGED*/

void meta_subset_columns(meta *src, meta *dest, int *ids, int nids) {

int i = 0;

  if (!(*dest).names && (*src).names)
    (*dest).names = Calloc1D((*src).ncols, sizeof(char *));
  if (!(*dest).flag && (*src).flag)
    (*dest).flag = Calloc1D((*src).ncols, sizeof(flags));

  for (i = 0; i < nids; i++) {

    if ((*src).names)
      (*dest).names[i] = (*src).names[ids[i]];
    if ((*src).flag)
      (*dest).flag[i] = (*src).flag[ids[i]];

  }/*FOR*/

  (*dest).nobs = (*src).nobs;
  (*dest).ncols = nids;

}/*META_SUBSET_COLUMNS*/

void FreeMETA(meta *m) {

  Free1D((*m).flag);
  Free1D((*m).names);

}/*FREEMETA*/

/* -------------------- discrete data ------------------------------------ */

ddata ddata_from_SEXP(SEXP df, int offset) {

int i = 0, ncols = length(df);
SEXP col;
ddata dt = empty_ddata(length(VECTOR_ELT(df, 0)), ncols + offset);

  for (i = 0; i < ncols; i++) {

    col = VECTOR_ELT(df, i);
    dt.col[i + offset] = INTEGER(col);
    dt.nlvl[i + offset] = NLEVELS(col);

  }/*FOR*/

  return dt;

}/*DDATA_FROM_SEXP*/

ddata empty_ddata(int nobs, int ncols) {

ddata dt = { 0 };

  dt.m.nobs = nobs;
  dt.m.ncols = ncols;
  dt.col = Calloc1D(ncols, sizeof(int *));
  dt.nlvl = Calloc1D(ncols, sizeof(int));

  return dt;

}/*EMPTY_DDATA*/

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
void FreeDDT(ddata dt, int free_data) {

  /* free the columns holding the data and the pointers, or free just the
   * pointers. */
  if (free_data)
    Free2D(dt.col, dt.m.ncols);
  else
    Free1D(dt.col);

  Free1D(dt.nlvl);
  FreeMETA(&(dt.m));

}/*FREEDDT*/

/* -------------------- Gaussian data ------------------------------------ */

/* create a data table from Gaussian data. */
gdata gdata_from_SEXP(SEXP df, int offset) {

int i = 0, ncols = length(df);
gdata dt = empty_gdata(length(VECTOR_ELT(df, 0)), ncols + offset);

  for (i = 0; i < ncols; i++)
    dt.col[i + offset] = REAL(VECTOR_ELT(df, i));

  return dt;

}/*GDATA_FROM_SEXP*/

gdata empty_gdata(int nobs, int ncols) {

gdata dt = { 0 };

  dt.m.nobs = nobs;
  dt.m.ncols = ncols;
  dt.col = Calloc1D(ncols, sizeof(double *));

  return dt;

}/*EMPTY_GDATA*/

/* cache the means of the continuous columns. */
void gdata_cache_means(gdata *dt, int offset) {

  (*dt).mean = Calloc1D((*dt).m.ncols, sizeof(double));
  c_meanvec((*dt).col, (*dt).mean, (*dt).m.nobs, (*dt).m.ncols, offset);

}/*CACHE_MEANS*/

void print_gdata(gdata dt) {

int i = 0;

  Rprintf("cgdata: %dx%d\n", dt.m.nobs, dt.m.ncols);

  for (i = 0; i < dt.m.ncols; i++) {

    Rprintf("%10s", dt.m.names[i]);
    Rprintf(" [%s%s]",
      (dt.m.flag[i].fixed ? "F" : " "),
      (dt.m.flag[i].drop ? "D" : " "));
    Rprintf("@%p", (void *)dt.col[i]);
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

}/*DDATA_SUBSET_COLUMNS*/

/* free a Gaussian data table. */
void FreeGDT(gdata dt, int free_data) {

  /* free the columns holding the data and the pointers, or free just the
   * pointers. */
  if (free_data)
    Free2D(dt.col, dt.m.ncols);
  else
    Free1D(dt.col);

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

  /* initialize the object. */
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

cgdata empty_cgdata(int nobs, int dcols, int gcols) {

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

  return dt;

}/*EMPTY_CGDATA*/

void print_cgdata(cgdata dt) {

int i = 0;

  Rprintf("cgdata: %dx%d\n", dt.m.nobs, dt.m.ncols);

  for (i = 0; i < dt.m.ncols; i++) {

    Rprintf("%10s", dt.m.names[i]);
    Rprintf(" [%s%s%s%s]",
      (dt.m.flag[i].discrete ? "D" : " "),
      (dt.m.flag[i].gaussian ? "G" : " "),
      (dt.m.flag[i].fixed ? "F" : " "),
      (dt.m.flag[i].drop ? "D" : " "));
    Rprintf("@%p", dt.m.flag[i].discrete ?
                     (void *)dt.dcol[dt.map[i]] : (void *)dt.gcol[dt.map[i]]);

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

void FreeCGDT(cgdata dt, int free_data) {

  /* free the columns holding the data and the pointers, or free just the
   * pointers. */
  if (free_data) {

    Free2D(dt.gcol, dt.ngcols);
    Free2D(dt.dcol, dt.ndcols);

  }/*THEN*/
  else {

    Free1D(dt.gcol);
    Free1D(dt.dcol);

  }/*ELSE*/

  Free1D(dt.nlvl);
  Free1D(dt.map);
  FreeMETA(&(dt.m));

}/*FREECGDT*/


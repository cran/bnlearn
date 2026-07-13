#ifndef DATA_TABLE_HEADER
#define DATA_TABLE_HEADER

/* flags for the columns (variables) of a data table, stored in the metadata. */
typedef struct {

  unsigned char own        : 1;  /* this column belongs with the struct. */
  unsigned char discrete   : 1;  /* this column is discrete data. */
  unsigned char continuous : 1;  /* this column is continuous data. */
  unsigned char complete   : 1;  /* this column contains no missing data points. */
  unsigned char fixed      : 1;  /* this column should never be (re)moved. */
  unsigned char drop       : 1;  /* this column is to be removed. */
  unsigned char padding    : 2;  /* pad to 1 byte. */

} flags;

/* columns (variables) metadata common to all types of data tables. */
typedef struct {

  int nobs;             /* number of observations. */
  int ncols;            /* number of columns. */
  const char **names;   /* column names (optional). */
  flags *flag;

} meta;

/* tabular data structure for data frames and their metadata. */
typedef struct {

  meta m;               /* metadata. */
  int **dcol;           /* pointers to the discrete columns. */
  double **ccol;        /* pointers to the continuous columns. */
  int *nlvl;            /* number of levels of the discrete columns. */
  double *mean;         /* means of the continuous columns (optional). */
  int ndcols;           /* number of discrete columns. */
  int nccols;           /* number of continuous columns. */
  int *map;             /* mapping between the original column position and the
                         * positions of the columns in the data structure. */

} tabular;

void meta_init_flags(meta *m, int offset, SEXP complete, SEXP fixed);
void meta_copy_names(meta *m, int offset, SEXP df);
void meta_copy(meta *src, meta *dest);
void FreeMETA(meta *m);

tabular tabular_from_SEXP(SEXP df, int doffset, int coffset);
tabular empty_tabular(int nobs, int dcols, int ccols);
tabular new_tabular(int nobs, int dcols, int ccols);
void tabular_cache_means(tabular *dt, int offset);
void print_tabular(tabular dt);
void tabular_drop_flagged(tabular *dt, tabular *copy);
void tabular_subset_columns(tabular *dt, tabular *copy, int *ids, int nids);
void tabular_incomplete_cases_range(tabular *dt, bool *indicator, int dcol_start,
    int dcol_end, int ccol_start, int ccol_end);
void tabular_incomplete_cases(tabular *dt, bool *indicator, int doffset,
    int coffset);
void tabular_subsample_by_logical(tabular *dt, tabular *copy, bool *indicators,
    int doffset, int coffset);
void FreeTAB(tabular dt);

#endif

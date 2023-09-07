#ifndef CONTINGENCY_TABLES_HEADER
#define CONTINGENCY_TABLES_HEADER

/* one-dimensional contingency table. */
typedef struct {

  int llx;     /* first (and only) dimension. */
  int nobs;    /* total count over all cells. */
  int *n;      /* contingency table. */

} counts1d;

/* two-dimensional contingency table. */
typedef struct{

  int llx;      /* first dimension. */
  int lly;      /* second dimension */
  int nobs;     /* total count over all cells. */
  int **n;      /* contingency table. */
  int *ni;      /* marginal counts for the first dimension. */
  int *nj;      /* marginal counts for the second dimension. */

} counts2d;

/* three-dimensional contingency table, as an array of two-dimensional tables
 * spanning the third dimension. */
typedef struct {

  int llx;      /* first dimension. */
  int lly;      /* second dimension */
  int llz;      /* third dimension. */
  int nobs;     /* total count over all cells. */
  int ***n;     /* contingency table. */
  int **ni;     /* marginal counts for the first dimension. */
  int **nj;     /* marginal counts for the second dimension. */
  int *nk;      /* marginal counts for the third dimension. */

} counts3d;

counts1d new_1d_table(int llx);
counts2d new_2d_table(int llx, int lly, bool margins);
counts3d new_3d_table(int llx, int lly, int llz);

void fill_1d_table(int *xx, counts1d *table, int num);
void fill_2d_table(int *xx, int *yy, counts2d *table, int num);
void fill_3d_table(int *xx, int *yy, int *zz, counts3d *table, int num);

void refill_1d_table(int *xx, counts1d *table, int num);
void refill_2d_table(int *xx, int *yy, counts2d *table, int num);
void refill_3d_table(int *xx, int *yy, int *zz, counts3d *table, int num);

void resize_1d_table(int llx, counts1d *table);
void resize_2d_table(int llx, int lly, counts2d *table);
void resize_3d_table(int llx, int lly, int llz, counts3d *table);

void print_1d_table(counts1d table);
void print_2d_table(counts2d table);
void print_3d_table(counts3d table);

void Free1DTAB(counts1d table);
void Free2DTAB(counts2d table);
void Free3DTAB(counts3d table);

void rcounts2d(counts2d table, double *fact, int *workspace);
void rcounts3d(counts3d table, double *fact, int *workspace);

#endif

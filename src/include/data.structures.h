#ifndef DATA_STRUCTURES_HEADER
#define DATA_STRUCTURES_HEADER

/* flags for the columns (variables) of a data table, stored in the metadata. */
typedef struct {

  unsigned int discrete : 1;  /* this column is discrete data. */
  unsigned int gaussian : 1;  /* this column is Gaussian data. */
  unsigned int complete : 1;  /* this column contains no missing data points. */
  unsigned int fixed    : 1;  /* this column should never be (re)moved. */
  unsigned int drop     : 1;  /* this column is to be removed. */
  unsigned int padding  : 3;  /* pad to 1 byte. */

} flags;

/* columns (variables) metadata common to all types of data tables. */
typedef struct {

  int nobs;             /* number of observations. */
  int ncols;            /* number of columns. */
  const char **names;   /* column names (optional). */
  flags *flag;

} meta;

/* data table for discrete data. */
typedef struct {

  meta m;               /* metadata. */
  int **col;            /* pointers to the discrete columns. */
  int *nlvl;            /* number of levels of the discrete columns. */

} ddata;

/* data table for gaussian data. */
typedef struct {

  meta m;               /* metadata. */
  double **col;         /* pointers to the continuous columns. */
  double *mean;         /* means of the continuous columns (optional). */

} gdata;

/* data table for conditional gaussian data. */
typedef struct {

  meta m;               /* metadata. */
  int **dcol;           /* pointers to the discrete columns. */
  double **gcol;        /* pointers to the continuous columns. */
  int *nlvl;            /* number of levels of the discrete columns. */
  int ndcols;           /* number of discrete columns. */
  int ngcols;           /* number of continuous columns. */
  int *map;             /* mapping between the original column position and the
                         * positions of the columns in the data structure. */

} cgdata;

/* covariance matrix, with additional fields to carry around its own SVD
 * decomposition and dimension. */
typedef struct {

  int dim;     /* dimension of the covariance matrix. */
  double *mat; /* pointer to the matrix. */
  double *u;   /* SVD decomposition, U left matrix, for c_svd(). */
  double *d;   /* SVD decomposition, vector for U's diagonal, for c_svd(). */
  double *vt;  /* SVD decomposition, V^t right matrix, for c_svd(). */

} covariance;

#endif

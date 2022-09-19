#ifndef COVARIANCE_MATRIX_HEADER
#define COVARIANCE_MATRIX_HEADER

/* covariance matrix, with additional fields to carry around its own SVD
 * decomposition and dimension. */
typedef struct {

  int dim;     /* dimension of the covariance matrix. */
  double *mat; /* pointer to the matrix. */
  double *u;   /* SVD decomposition, U left matrix, for c_svd(). */
  double *d;   /* SVD decomposition, vector for U's diagonal, for c_svd(). */
  double *vt;  /* SVD decomposition, V^t right matrix, for c_svd(). */

} covariance;

covariance new_covariance(int dim, bool decomp);
void print_covariance(covariance cov);
void copy_covariance(covariance *src, covariance *copy);
void covariance_drop_variable(covariance *full, covariance *sub, int to_drop);
void FreeCOV(covariance cov);

void c_covmat(double **data, double *mean, int nrow, int ncol, covariance cov,
    int first);
void c_update_covmat(double **data, double *mean, int update, int nrow,
    int ncol, double *mat);
void c_covmat_with_missing(double **data, int nrow, int ncol,
    bool *missing_partial, bool *missing_all, double *mean, double *mat,
    int *ncomplete);

#endif

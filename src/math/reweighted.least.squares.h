#ifndef REWEIGHTED_LEAST_SQUARES_HEADER
#define REWEIGHTED_LEAST_SQUARES_HEADER

/* GLM family for IRLS. */
typedef enum {
  GLM_NEGBIN = 0,        /* negative binomial (NB2 parameterisation). */
  GLM_HYPERPOISSON = 1   /* hyper-Poisson. */
} glm_family_e;

/* scratch buffers for one IRLS fit, allocated once by the caller and reused
 * across calls (w, z, zs: length n; beta_new: length p; Xs: length n * p). */
typedef struct {
  double *w, *z, *beta_new, *Xs, *zs;
} irls_scratch;

/* link functions. */
double irls_inv_logit(double eta);
double irls_linpred(const double *X, int n, int p, const double *coef, int i);

int c_weighted_least_squares(const double *X, int n, int p, const double *w,
    const double *z, double *beta, double *Xs, double *zs);

/* weighted logistic regression (binomial family, logit link) by IRLS. */
int c_irls_logistic(const double *X, int n, int p, const double *y,
    const double *prior_w, double *beta, int max_iter, double tol,
    irls_scratch *ws);

/* weighted hyper-Poisson or negative binomial regression (log link) by IRLS. */
int c_irls_count(glm_family_e family, const double *X, int n, int p,
    const double *y, const double *prior_w, double *beta, double *dispersion,
    int max_iter, double tol, irls_scratch *ws);

#endif

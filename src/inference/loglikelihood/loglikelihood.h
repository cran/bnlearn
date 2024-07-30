#ifndef LOGLIKELIHOOD_FUNCTIONS_HEADER
#define LOGLIKELIHOOD_FUNCTIONS_HEADER

bool check_locally_incomplete_data(fitted_bn bn, meta m, bool debugging);

void bysample_discrete_loglikelihood(fitted_bn bn, ddata dt, double *loglik,
    bool debugging);
void bysample_gaussian_loglikelihood(fitted_bn bn, gdata dt, double *loglik,
    bool robust, bool debugging);
void bysample_clgaussian_loglikelihood(fitted_bn bn, cgdata dt, double *loglik,
    bool robust, bool debugging);

double data_discrete_loglikelihood(fitted_bn bn, ddata dt, bool propagate,
    bool loss, bool debugging);
double data_gaussian_loglikelihood(fitted_bn bn, gdata dt, double *scratch,
    bool loss, bool propagate, bool debugging);
double data_clgaussian_loglikelihood(fitted_bn bn, cgdata dt, double *scratch,
    bool loss, bool propagate, bool debugging);

#endif


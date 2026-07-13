#ifndef PREDICTIONS_HEADER
#define PREDICTIONS_HEADER

#include "../core/data.table.h"
#include "../fitted/fitted.h"
#include "../include/sampling.h"

void c_naivepred(tabular dt, fitted_bn bn, int tr_id, int tr_nlvl, double *pr,
    int *map, int *res, double *pt, bool include_prob, bool debugging);
void c_mappred(tabular dt, fitted_bn bn, int node_id, bool target_discrete,
    int nlvls, int nev, bool *ev_discrete, void **varptrs, void **evptrs,
    double *wgt, int nsims, long double *lvls_counts, void *pred, void *res,
    double *pt, bool include_prob, int *drop, bool debugging, SEXP fitted,
    SEXP cpdist, SEXP from, fixed_node *fixed);

/* SEXP-free kernels for parents-based predictions (predict/parents.c). */
void c_gaussian_predict(ldist ld, double **ccols, int ncont, int nobs,
    double *res, bool debugging);
void c_discrete_predict(ldist ld, int *configs, int nconfigs, int nobs,
    int *res, double *pt, bool include_prob, bool debugging);
void c_cgaussian_predict(ldist ld, double **ccols, int ncont, int *configs,
    int nobs, double *res, bool debugging);
void c_zero_inflated_predict(ldist ld, fitted_node_e type, tabular parents, int nobs,
    double *res, bool debugging);

#endif

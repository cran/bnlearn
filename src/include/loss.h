
/* from loss.c */
double c_entropy_loss(SEXP fitted, SEXP orig_data, int ndata, int by,
    double *res_sample, SEXP keep, bool allow_singular, bool warn,
    bool debugging);


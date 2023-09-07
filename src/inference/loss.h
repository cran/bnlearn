
/* from loss.c */
double c_entropy_loss(SEXP fitted, SEXP orig_data, int ndata, bool by,
    double *res_sample, double *effective, SEXP keep, bool allow_singular,
    bool warn, bool debugging);


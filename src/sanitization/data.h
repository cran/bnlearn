#ifndef DATA_SANITIZATION_HEADER
#define DATA_SANITIZATION_HEADER

SEXP data_type(SEXP data);
SEXP count_observed_values(SEXP data);
void add_simulation_metadata(SEXP result, int num);

#endif

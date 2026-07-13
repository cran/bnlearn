#ifndef HYPERPOISSON_HEADER
#define HYPERPOISSON_HEADER

double dhypois(double x, double intensity, double dispersion, bool give_log);
double rhypois(double intensity, double dispersion);
void hypois_moments(double intensity, double dispersion, double *mu, double *var);
double hypois_mean(double intensity, double dispersion);

#endif

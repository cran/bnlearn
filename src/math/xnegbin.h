#ifndef XNEGATIVEBINOMIAL_HEADER
#define XNEGATIVEBINOMIAL_HEADER

double dxnegbin(double x, double p, double fails, bool give_log);
double negbin_mean(double p, double fails);
double rxnegbin(double p, double fails);

#endif

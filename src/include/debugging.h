
/* instrumentation to debug PROTECT()/UNPROTECT() calls. */
void PROTECT_DEBUG(SEXP s, const char *fun, const char *file, int line);
void UNPROTECT_DEBUG(int n, const char *fun, const char *file, int line);
#ifdef DEBUG
#undef PROTECT
#define PROTECT(s) PROTECT_DEBUG(s, __func__, __FILE__, __LINE__)
#undef UNPROTECT
#define UNPROTECT(n) UNPROTECT_DEBUG(n, __func__, __FILE__, __LINE__)
#endif



/* the global test counter, now living in C-land. */
extern double test_counter;

/* symbols used to set object attributes. */
extern SEXP BN_ModelstringSymbol;
extern SEXP BN_NodesSymbol;
extern SEXP BN_ProbSymbol;
extern SEXP BN_MethodSymbol;
extern SEXP BN_WeightsSymbol;
extern SEXP BN_DsepsetSymbol;
extern SEXP BN_MetaDataSymbol;
extern SEXP BN_NobsSymbol;
extern SEXP BN_DfSymbol;

/* shortcuts to boolean SEXPs. */
extern SEXP TRUESEXP, FALSESEXP;

/* numerical constants */
#define MACHINE_TOL sqrt(DBL_EPSILON)


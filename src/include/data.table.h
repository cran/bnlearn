#include "data.structures.h"

void meta_init_flags(meta *m, int offset, SEXP complete, SEXP fixed);
void meta_copy_names(meta *m, int offset, SEXP df);

ddata ddata_from_SEXP(SEXP df, int offset);
ddata empty_ddata(int nobs, int ncols);
void ddata_drop_flagged(ddata *dt, ddata *copy);
void ddata_subset_columns(ddata *dt, ddata *copy, int *ids, int nids);
void FreeDDT(ddata dt, int free_data);

gdata gdata_from_SEXP(SEXP df, int offset);
gdata empty_gdata(int nobs, int ncols);
void gdata_cache_means(gdata *dt, int offset);
void gdata_move_column(gdata *dt, gdata *copy, int i, int j);
void print_gdata(gdata dt);
void gdata_drop_flagged(gdata *dt, gdata *copy);
void gdata_subset_columns(gdata *dt, gdata *copy, int *ids, int nids);
void FreeGDT(gdata dt, int free_data);

cgdata cgdata_from_SEXP(SEXP df, int doffset, int goffset);
cgdata empty_cgdata(int nobs, int dcols, int gcols);
void print_cgdata(cgdata dt);
void cgdata_drop_flagged(cgdata *dt, cgdata *copy);
void cgdata_subset_columns(cgdata *dt, cgdata *copy, int *ids, int nids);
void FreeCGDT(cgdata, int free_data);


/* memory allocation functions */
void *R_calloc(size_t num, size_t size);
#define alloc1dcont(x) ((int *) R_calloc(x, sizeof(int)))
int **alloc2dcont(int length, int width);
int ***alloc3dcont(int length, int width, int depth);
#define allocstatus(x) ((short int *) R_calloc(x, sizeof(short int)))
#define alloc1dreal(x) ((double *) R_calloc(x, sizeof(double)))
double **alloc2dreal(int length, int width);
double ***alloc3dreal(int length, int width, int depth);
#define alloc1dpointer(x) ((void **) R_calloc(x, sizeof(void *)))
#define alloc1dstring(x) ((char **) alloc1dpointer(x))
#define allocldouble(x) ((long double *) R_calloc(x, sizeof(long double)))


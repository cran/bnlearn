
/* memory allocation functions */
void *R_calloc(size_t num, size_t size);
void **R_calloc2(size_t num1, size_t num2, size_t size);
void ***R_calloc3(size_t num1, size_t num2, size_t num3, size_t size);
#define alloc1dcont(x) ((int *) R_calloc(x, sizeof(int)))
#define alloc2dcont(x, y) ((int **) R_calloc2(x, y, sizeof(int)))
#define alloc3dcont(x, y, z) ((int ***) R_calloc3(x, y, z, sizeof(int)))
#define allocstatus(x) ((short int *) R_calloc(x, sizeof(short int)))
#define alloc1dreal(x) ((double *) R_calloc(x, sizeof(double)))
#define alloc2dreal(x, y) ((double **) R_calloc2(x, y, sizeof(double)))
#define alloc3dreal(x, y, z) ((double ***) R_calloc3(x, y, z, sizeof(double)))
#define alloc1dpointer(x) ((void **) R_calloc(x, sizeof(void *)))
#define alloc1dstring(x) ((char **) alloc1dpointer(x))
#define allocldouble(x) ((long double *) R_calloc(x, sizeof(long double)))


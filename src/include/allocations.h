
/* memory allocation functions */
int *alloc1dcont(int length);
int **alloc2dcont(int length, int width);
int ***alloc3dcont(int length, int width, int depth);
short int *allocstatus(int length);
double *alloc1dreal(int length);
double **alloc2dreal(int length, int width);
double ***alloc3dreal(int length, int width, int depth);
void **alloc1dpointer (int length);
char **alloc1dstring (int length);
long double *allocldouble (int length);


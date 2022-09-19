#ifndef ALLOCATIONS_HEADER
#define ALLOCATIONS_HEADER

/* memory allocation. */
void *Calloc1D(size_t R, size_t size);
void **Calloc2D(size_t R, size_t C, size_t size);
void ***Calloc3D(size_t R, size_t C, size_t L, size_t size);
void *Realloc1D(void *p, size_t R, size_t size);
void BN_Free1D(void *p);
void BN_Free2D(void **p, size_t R);
void BN_Free3D(void ***p, size_t R, size_t C);

#define Free1D(p) \
  do { \
    BN_Free1D((void *)p); \
    p = NULL; \
  } while(0)
#define Free2D(p, R) \
  do { \
    BN_Free2D((void **)p, R); \
    p = NULL; \
  } while (0)
#define Free3D(p, R, C) \
  do { \
    BN_Free3D((void ***)p, R, C); \
    p = NULL; \
  } while (0)

#endif

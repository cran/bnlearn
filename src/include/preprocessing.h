
/* enum for discretization methods, to be matched from the label string passed
 * down from R. */
typedef enum {
  ENOMETHOD    =  0, /* error code, no such discretization method. */
  INTERVAL     =  1, /* interval discretization (marginal). */
  QUANTILE     =  2, /* quantile discretization (marginal). */
  HARTEMINK    =  3  /* Hartemink discretization (pairwise). */
} discretization_e;

discretization_e discretization_to_enum(const char *label);


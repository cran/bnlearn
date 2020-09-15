#include "data.structures.h"

/* from contingency.tables.c */
counts1d new_1d_table(int llx);
counts2d new_2d_table(int llx, int lly, bool margins);
counts3d new_3d_table(int llx, int lly, int llz);
void fill_1d_table(int *xx, counts1d *table, int num);
void fill_2d_table(int *xx, int *yy, counts2d *table, int num);
void fill_3d_table(int *xx, int *yy, int *zz, counts3d *table, int num);
void Free1DTAB(counts1d table);
void Free2DTAB(counts2d table);
void Free3DTAB(counts3d table);

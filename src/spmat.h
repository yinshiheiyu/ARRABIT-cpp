#include "mkl.h"

typedef struct {
    double *val;
    MKL_INT *indx;
    MKL_INT *pntr;
    int nnz;
    int col;
} SpMat;
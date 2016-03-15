#include "mkl.h"

typedef struct {
    double *val;
    MKL_INT *indx;
    MKL_INT *pntr;
    int nnz;
    int row;
    int col;
} SpMat;
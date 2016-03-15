#include<stdio.h>
#include "mkl.h"
#include "spmat.h"

void read_csc(const char *filename, SpMat* A){
    FILE *fp = fopen(filename, "rb");
    int nnz, col;
    if(fp == NULL){
        printf("Cannot open file %s. \n", filename);
        return;
    }
    fread(&nnz, sizeof(MKL_INT), 1, fp);
    A->nnz = nnz;
    A->val = (double*)malloc(nnz * sizeof(double));
    A->indx = (MKL_INT*)malloc(nnz * sizeof(MKL_INT));
    fread(A->val, sizeof(double), nnz, fp);
    fread(A->indx, sizeof(MKL_INT), nnz, fp);
    fread(&col, sizeof(MKL_INT), 1, fp);
    A->col = col;
    A->pntr = (MKL_INT*)malloc((col + 1) * sizeof(MKL_INT));
    fread(A->pntr, sizeof(MKL_INT), col + 1, fp);
    fclose(fp);
}
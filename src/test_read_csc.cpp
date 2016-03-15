#include<stdio.h>
#include "spmat.h"

void read_csc(const char*, SpMat*);

int main(){
    char *filename = "../matbin/csc/Ga3As3H12";
    SpMat A;
    read_csc(filename, &A);
    printf("nnz : %d\n", A.nnz);
    printf("val : %p\n", A.val);
    printf("indx : %p\n", A.indx);
    printf("col : %d\n", A.col);
    printf("pntr : %p\n", A.pntr);
    return 0;
}
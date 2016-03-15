#include "mkl.h"
#include <omp.h>

/* norm1 - computes the 1-norm of an m-by-n matrix
 *   Inputs :
 *         m, n - size of A
 *         A    - array, ColMajorLayout
 *   Return Value:
 *         (double) - the 1-norm of the matrix
 */
double norm1(double* A, const MKL_INT m, const MKL_INT n){
	double *sums = new double[n], max = 0;
	int i;
	CBLAS_INDEX indx;
#	pragma omp parallel for num_threads(Global::thread_num)\
	default(none) private(i) shared(A, sums)
	for(i = 0; i < n; i++){
		sums[i] = cblas_dasum(m, A + i * m, 1);
	}
	indx = cblas_idamax(n, sums, 1);
	max = sums[indx];
	delete []sums;
	return max;
}

/* normCol2 - normalize each column of an m-by-n matrix (2-norm)
 *   Inputs:
 *         A    -  array, represents the m-by-n matrix, ColMajorLayout
 *         m, n -  size of A, A is m-by-n
 */
void normCol2(double *A, const MKL_INT m, const MKL_INT n){
	int i;
	double local_res;
#	pragma omp parallel for num_threads(Global::thread_num) \
	default(none) private(i, local_res) shared(A)
	for(i = 0; i < n; i++){
		local_res = cblas_dnrm2(m, A + i * m, 1);
		if(local_res != 0){
			// cblas_dscal : x := ax
			cblas_dscal(m, 1 / local_res, A + i * m, 1);
		}
	}
}
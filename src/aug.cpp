#include "arrabit_inc.h"

/* doAug - do augmentation Y = [X, AX, ..., A^p X]
 *   Inputs :
 *         n, k                          - A is n-by-n; X is n-by-k
 *         val, indx, pointerB, pointerE - represents the sparse matrix A
 *         X                             - input matrix X
 *         nblks                         - number of blocks, nblks = p + 1
 *         ln                            - ???
 *         Z                             - pre-allocated matrix Z, Z is n-by-(nblks*k)
 */
void doAug(const MKL_INT n, const MKL_INT k, const double* val, const MKL_INT* indx, const MKL_INT* pointerB, const MKL_INT* pointerE, const double* X, const int nblks, double ln, double* Z){
	char matdescra[] = {'G', 'L', 'N', 'F'};
	char trans[] = {'N'};
	double alpha = 1.0, beta = -ln;
	double *Y = new double[n * k];
	mat_dcopy(n, k, X, Y);
	mat_dcopy(n, k, X, Z);
	for(int i = 1; i < nblks; i++){
		//Y = A * Y - ln * Y
		mkl_dcscmm(trans, &n, &k, &n, &alpha, matdescra, val, indx, pointerB, pointerE, Y, &n, &beta, Y, &n);
		//Z = [Z Y]
		mat_dcopy(n, k, Y, Z + i * n * k);
	}
	//normalize
	normCol2(Z, n, k * nblks);
	delete []Y;
}
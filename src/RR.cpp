#include "arrabit_inc.h"

/* RR - Alg.1 do RR projection [U, d] = RR(A, Z)
 *   Inputs :
 *         n, k                          - input Z is n-by-k
 *         val, indx, pointerB, pointerE - represents the sparse matrix A
 *         Z                             - array, size n-by-k
 *         U                             - pre-allocated array, size at least n-by-k, serves as the evc output
 *         d                             - pre-allocated array, size at least k, serves as the ev output
 *   Notes  :
 *         array Z will be changed after executing this function
 */
void RR(const MKL_INT n, const MKL_INT k, const double* val, const MKL_INT* indx, const MKL_INT* pointerB, const MKL_INT* pointerE, double* Z, double* U, double* d){
	char matdescra[] = {'G', 'L', 'N', 'F'};
	char trans[] = {'N'};
	double alpha = 1.0, beta = 0.0, tmp;
	int i;
	//perform QR factorization
	double tau = new double[k];
	LAPACKE_dgeqrf(LAPACK_COL_MAJOR, n, k, Z, n, tau);
	//extract Q[1:k]
	dorgqr(LAPACK_COL_MAJOR, n, k, k, Z, n, tau);
	double AZ = new double[n * k];
	double ZtAZ = new double[k * k];
	double V = new double[k * k];
	//AZ = A * Z
	mkl_dcscmm(trans, &n, &k, &n, &alpha, matdescra, val, indx, pointerB, pointerE, Z, &n, &beta, AZ, &n);
	//ZtAZ = Z' * AZ
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, k, k, n, 1.0, Z, n, AZ, n, 0.0, ZtAZ, k);
	//compute eigenvalue factorization
	eigqr(k, ZtAZ, V, d);
	//array d contains all eigenvalues in ascending order, thus convert d and V into descending order
#	pragma omp parallel for num_threads(OMP_THREAD_NUM) \
	default(none) private(i, tmp) shared(V, d, k, n)
	for(i = 0; i < k / 2;i++){
		tmp = d[k - i - 1]; d[k - i - 1] = d[i]; d[i] = tmp;
		//since V is stored in ColMajor, swap V(:,i) and V(:, k - 1 - i)
		cblas_dswap(n, V + i * n, 1, V + (k - 1 - i) * n, 1);
	}
	//U = Z * V
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, k, k, 1.0, Z, n, V, k, 0.0, U, n);
}
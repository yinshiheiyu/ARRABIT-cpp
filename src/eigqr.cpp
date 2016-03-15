#include "arrabit_inc.h"

/* eigqr - compute all eigenvalues and eigenvectors of symmetric matrix A
 *   Inputs :
 *         n - dimension of matrix A
 *         A - array, size n-by-n
 *         V - pre-allocated array, size at least n-by-n, stores all eigenvectors
 *         d - pre-allocated array, size at least n, stores all eigenvalues
 */
void eigqr(const MKL_INT n, double* A, double* V, double* d){
	double *e, *tau;
	e = new double[n];
	tau = new double[n];
	//reduce A to tridiagonal form
	LAPACKE_dsytrd(LAPACK_COL_MAJOR, 'U', n, A, n, d, e, tau);
	//compute all eigenvalues and eigenvectors via QR alg.
	LAPACKE_dsteqr(LAPACK_COL_MAJOR, 'I', n, d, e, V, n);
	//obtain real V, i.e. V <- QV
	LAPACKE_dormtr(LAPACK_COL_MAJOR, 'L', 'U', 'N', n, n, A, n, tau, V, n);
	delete []e;
	delete []tau;
}
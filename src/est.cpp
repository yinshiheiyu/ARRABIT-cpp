#include "mkl.h"
#include "spmat.h"

/* est_key_eigenvalues - estimate key eigenvalues l1, lkw, ln
 *   Input :
 *        n, kw                   - X is n-by-kw, A is n-by-n
 *        X                       - array, represents the matrix X
 *        val, indx, pntrB, pntrE - represents the matrix X
 *        which                   - integer 0~7, indicates which key eigenvalue is given, l1-lkw-ln
 *        maxit                   - max iterations during RR procedure
 *        l1, lkw, ln             - initial eigenvalues
 *   Output :
 *        l1, lkw, ln             - estimated key eigenvalues
 */
void est_key_eigenvalues(const MKL_INT n, const MKL_INT kw, double* val, MKL_INT* indx, MKL_INT* pntrB, MKL_INT* pntrE, double* X, int which, int maxit, double* l1, double* lkw, double* ln){
	int ierr, i;
	char sa[] = {'S', 'A'}, be[] = {'B', 'E'}, trans[] = {'N'}, matdescra[] = {'G', 'L', 'N', 'F'};
	double *Y, w[2], alpha = 1.0, ln_neg, tau, *T;
	//estimate ln
	if(which & 0x01 == 0){ //ln is not given
		ierr = eigs_sps(n, val, indx, pntrB, pntrE, sa, 1, ln);
	}
	ln_neg = -ln[0];
	//return if lkw and l1 are both given
	if(which & 0x02 && which & 0x04) return;
	Y = new double[n * kw];
	//power method (roughly)
	for(i = 0; i < maxit; i++){
		//Y := A * Y - ln * Y
		mkl_dcscmm(trans, &n, &kw, &n, &alpha, matdescra, val, indx, pntrB, pntrE, Y, &n, &ln_neg, Y, &n);
	}
	//RR procedure
	tau = new double[kw];
	LAPACKE_dgeqrf(LAPACK_COL_MAJOR, n, kw, Y, n, tau);
	//extract Y[1:kw]
	dorgqr(LAPACK_COL_MAJOR, n, kw, kw, Y, n, tau);
	//X = A * Y - ln * Y
	mat_dcopy(n, kw, Y, X);
	mkl_dcscmm(trans, &n, &kw, &n, &alpha, matdescra, val, indx, pntrB, pntrE, X, &n, &ln_neg, X, &n);
	//YtAY = Y' * X
	T = new double[kw * kw];
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, kw, kw, n, 1.0, Y, n, X, n, 0.0, T, kw);
	//estimations from T
	if(which & 0x04 == 0){ //l1 is not given
		eigs_ges(kw, T, kw, be, 2, w);
	}
	else{ //estimate lkw only
		eigs_ges(kw, T, kw, sa, 1, w);
	}
	*l1 = w[1] + *ln;
	*lkw = w[0] + *ln;
	delete []Y;
	delete []tau;
	delete []T;
}
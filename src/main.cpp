#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include "arrabit_inc.h"

int OMP_THREAD_NUM = 1;
int main(int argc, char **argv){
	return 0;
}
/* arrabit - the main alg.
 *   Inputs:
 *         n                     - size of A
 *         k                     - numbers of Ritz pairs needed to be computed
 *         A, indx, pntrB, pntrE - represents the sparse matrix A
 *         sigma                 - 'L' for Largest Algebraic; 'S' for Smallest Algebraic
 *         opts                  - optional parameters, using default settings when NULL
 *   Outputs:
 *         U                     - computed k extreme eigenvectors
 *         D                     - computed k extreme eigenvalues
 *   Return Value:
 *         (int)                 - 0  : process exits normally
 *                               - -1 : process exits with errors
 */
int arrabit(const MKL_INT n, const MKL_INT k, const double* A, const MKL_INT* indx, const MKL_INT* pntrB, const MKL_INT* pntrE, const char sigma, Opts* opts, double* U, double *D){
	int i;
	//begin parameters initialize
	int maxit0, maxit1, maxit2, cheby;
	char inner_solver;
	double tol, tol_t, tol_d, guard;
	double tiny, tiny1, tiny2;
	if(opts == NULL){
		maxit0 = 5; maxit1 = 10; maxit2 = 30; cheby = 1;
		inner_solver = 'M';
		guard = 0.1;
		tol = tol_t = 1E-6;
		tiny = 1E-14, tiny1 = 1E-13; tiny2 = 1E-12;
	}
	else{
	}
	if(tol_t > 1E-10) tol_d = tiny2;
	else if(tol_t > 1E-12) tol_d = tiny1;
	else tol_d = tiny > tol_t * tol_t ? tiny : tol_t * tol_t;
	//end parameters initialize
	MKL_INT kw = k + ceil(k * guard);
	MKL_INT deg;
	double a, b;
	double *X = new double[n * kw];
	double *X_best = new double[n * kw];
	double *Y = new double[n * kw];
	double *W = new double[kw * kw];
	double *ev = new double[kw];
	double *ev_best = new double[kw];
	double *res = new double[k];
	double rc, rcp, rc2, anorm;
	double tol_t, tol_d, maxres, maxres_prev = 0, bestres;
	int count_no_decay = 0, num_converged = 0;
	lapack_int *ipiv = new lapack_int[kw];
	normRnd(n, kw, X, 0, 1);
	for(int i1 = 0; i1 < maxit2; i1++){ //outer loop
		for(int i2 = 0; i2 < maxit1; i2++){ //level-1 loop
			rc = 10;
			for(int i3 = 0; i3 < maxit0; i3++){ //level-0 loop
				switch(inner_solver){
				case 'M': //MPM
					if(cheby) polyCheby(n, kw, A, indx, pntrB, pntrE, X, a, b, deg, Y);
					else polyFunc(); //??
					mat_dcopy(n, kw, Y, X);
					normCol2(X, n, kw);
					//proj2Ker??
					break;
				case 'G': //GN
					break;
				default:
					printf("Invalid inner-solver configuration...Ending process...\n");
					return -1;
				}
			} // end level-0 loop
			//W = X'X
			cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, kw, kw, n, 1.0, X, n, X, n, 0.0, W, kw);
			//compute rcond(W)
			anorm = norm1(W, kw, kw);
			LAPACKE_dgetrf(LAPACK_COL_MAJOR, kw, kw, W, kw, ipiv); //LU factorization
			LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', (lapack_int)kw, W, (lapack_int)kw, anorm, &rc);
			rcp = rc; rc2 = rc/rcp;
			if(rc < tol || rc2 > 0.99) break;
			//
		} // end level-1 loop
		//Extract leading Ritz pairs
		//compute residuals
		residual(n, k, A, indx, pntrB, pntrE, X, ev, res);
		//record the best residuals
		if(i1 == 0 || maxres < bestres){
			bestres = maxres;
			mat_dcopy(n, kw, X, X_best);
			cblas_dcopy(kw, ev, 1, ev_best, 1);
		}
		//check outer loop stop rule
		if(i1 > 0) maxres_prev = maxres;
		maxres = res[cblas_idamax(k, res, 1)];
		if(i1 > 0 && maxres > maxres_prev) count_no_decay += 1;
		else count_no_decay = 0;
		num_converged = 0;
#		pragma omp parallel for num_threads(OMP_THREAD_NUM)\
		default(none) private(i) shared(res) reduction(+:num_converged)
		for(i = 0; i < k; i++){
			if(res[i] < 0.1 * tol) num_converged += 1;
		}
		if(i1 == maxit2 - 1 || count_no_decay > 2){
			printf("The result may be incorrect.\n");
			mat_dcopy(n, k, X_best, U);
			cblas_dcopy(k, ev_bset, 1, D, 1);
			return -1;
		}
		if(maxres < tol || maxres < (1 + 9.0 * num_converged / k) * tol){
			mat_dcopy(n, k, X_best, U);
			cblas_dcopy(k, ev_best, 1, D, 1);
			return 0;
		}
		//continuation
		if(maxres < tol_t || maxres < (1 + 9.0 * num_converged / k) * tol_t){
			tol_t *= 0.01;
			if(tol_t < tol) tol_t = tol;
			b = ev[kw - 1];
		}
		//update augmentation blocks num
		if(ev[kw - 1] / ev[k - 1] > 0.95 && maxres / maxres_p > 0.1){
			nblks += 1;
			if(nblks > maxblks) nblks = maxblks;
		}			
		//deflation??
		if(tol_t > 1E-10) tol_d = tiny2;
		else if(tol_t > 1E-12) tol_d = tiny1;
		else tol_d = tiny > tol_t * tol_t ? tiny : tol_t * tol_t;
		//update polynomial degree
	} // end outer loop
}

/* mat_dcopy - copy matrix data from dest to target
 *   Inputs:
 *         m, n   -  size of dest, m must be the leading dimention of dest
 *         src    -  array, represents the input matrix, ColMajorLayout
 *         dest   -  pre-allocated array, size at least m-by-n
 */
void mat_dcopy(const MKL_INT m, const MKL_INT n, const double *dest, double *target){
	cblas_dcopy(m * n, src, 1, dest, 1);
}

/* proj2Qker - projection to [Qt Qc], i.e, Y = Y - Q(Q'Y)
 *   Inputs :
 *         n          - number of rows of Y, Qt, Qc
 *         k          - number of columns of Y
 *         k_Qt, k_Qc - number of columns of Qt and Qc respectively
 *         Y          - array, size at least n-by-k
 *         Qt, Qc     - array, size at least n-by-k_Qt, n-by-k_Qc respectively
 *   Notes  :
 *         array Y will be changed after executing this function
 */
void proj2Qker(const MKL_INT n, const MKL_INT k, const MKL_INT k_Qt, const MKL_INT k_Qc, double* Y, const double* Qt, const double* Qc);
	double *Q, *QtY;
	if(k_Qt + k_Qc == 0) return;
	Q = new double[n * (k_Qt + k_Qc)];
	if(k_Qt > 0) mat_dcopy(n, k_Qt, Qt, Q);
	if(k_Qc > 0) mat_dcopy(n, k_Qc, Qc, Q + n * k_Qt);
	//QtY : Q'Y
	QtY = new double[(k_Qt + k_Qc) * k];
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, k_Qt + k_Qc, k, n, 1.0, Q, n, Y, n, 0.0, QtY, k_Qt + K_Qc);
	//Y = Y - Q*QtY
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, k, k_Qt + k_Qc, -1.0, Q, n, QtY, k_Qt + k_Qc, 1.0, Y, n);
	delete []Q;
	delete []QtY;
}

#include "mkl.h"
#include <omp.h>

extern int OMP_THREAD_NUM;

/* polyAX - Alg.7 Polynomial function Y=p(A)X, Horner alg.
 *   Inputs :
 *         n, k                          - A is n-by-n; X is n-by-k
 *         val, indx, pointerB, pointerE - represents the sparse matrix A
 *         X                             - dense matrix X
 *         a, b                          - interval [a, b]
 *         p                             - array, represents the polynomial, p[0]~p[deg]
 *         deg                           - the degree of the polynomial
 *         Y                             - pre-allocated matrix Y
 *   Outputs:
 *         Y - array, m-by-n matrix
 */
void polyAX(const MKL_INT n, const MKL_INT k, const double* val, const MKL_INT* indx, const MKL_INT* pointerB, const MKL_INT* pointerE, const double *X, double a, double b, const double* p, const MKL_INT deg, double *Y){
	double c0 = (a + b) / (a - b), c1 = 2 / (b - a);
	char matdescra[] = {'G', 'L', 'N', 'F'};
	char trans[] = {'N'}
	int i;
	cblas_dcopy(k * n, X, 1, Y, 1);
	cblas_dscal(k * n, p[0], Y, 1);
	for(i = 1; i <= deg; i++){
		//Y = c1 * AY + c0 * Y
		mkl_dcscmm(trans, &n, &k, &n, &c1, matdescra, val, indx, pointerB, pointerE, Y, &n, &c0, Y, &n);
		//Y = Y + p[i] * X
		cblas_daxpy(k * n, p[i], X, 1, Y, 1);
	}
}

/* polyChebyAX - Alg.7 Polynomial function Y=p(A)X, Horner alg, using Chebyshev polynomials of degree deg
 *   Inputs :
 *         n, k                          - A is n-by-n; X is n-by-k
 *         val, indx, pointerB, pointerE - represents the sparse matrix A
 *         X                             - dense matrix X
 *         a, b                          - interval [a, b]
 *         deg                           - the degree of the polynomial
 *         Y                             - pre-allocated matrix Y
 *   Outputs:
 *         Y - array, m-by-n matrix
 */
void polyChebyAX(const MKL_INT n, const MKL_INT k, const double* val, const MKL_INT* indx, const MKL_INT* pntrB, const MKL_INT* pntrE, const double* X, double a, double b, const MKL_INT deg, double *Y){
	double beta = 0.0, alpha = 1.0, c0 = (a + b) / (a - b), c1 = 2 / (b - a), c02 = 2 * c0, c12 = 2 * c1;
	double *T1, *T2;
	char matdescra[] = {'G', 'L', 'N', 'F'};
	char trans[] = {'N'}
	int i;
	cblas_dcopy(k * n, X, 1, Y, 1);
	if(deg == 0) return;
	if(deg == 1){
		//Y = AY
		mkl_dcscmm(trans, &n, &k, &n, &alpha, matdescra, val, indx, pointerB, pointerE, Y, &n, &beta, Y, &n);
		return;
	}
	T1 = new double[n * k];
	T2 = new double[n * k];
	//T1 = X
	cblas_dcopy(k * n, X, 1, T1, 1);
	//T2 = c0 * X + c1 * A * X
	cblas_dcopy(k * n, X, 1, T2, 1);
	mkl_dcscmm(trans, &n, &k, &n, &c1, matdescra, val, indx, pointerB, pointerE, T2, &n, &c0, T2, &n);
	for(i = 2; i <= deg; i++ ){
		//Y = T2
		cblas_dcopy(k * n, T2, 1, Y, 1);
		//Y = 2 * c0 * Y + 2 * c1 * A * Y
		mkl_dcscmm(trans, &n, &k, &n, &c12, matdescra, val, indx, pointerB, pointerE, Y, &n, &c02, Y, &n);
		//Y = Y - T1
		cblas_daxpy(k * n, -1, T1, 1, Y, 1);
		if(i < deg){
			//T1 = T2
			cblas_dcopy(k * n, T2, 1, T1, 1);
			//T2 = Y
			cblas_dcopy(k * n, Y, 1, T2, 1);
		}
	}
	delete []T1;
	delete []T2;
}

/* ployValue - compute polynomial p(X), Horner Alg, p() is applied to each element of X
 *   Inputs :
 *         n    - the size of X
 *         X    - array, size at least n
 *         a, b - interval [a, b]
 *         pc   - array, polynomial coefficients, p[0]~p[deg]
 *         deg  - the degree of the polynomial
 *   Outputs:
 *         Y    - preallocated array, size at least n
 */
void polyValue(const MKL_INT n, const double* X, double a, double b, const double* pc, const MKL_INT deg, double* Y){
	double c0 = (a + b) / (a - b), c1 = 2 / (b - a);
	double *Bx = new double[n];
	double *ones = new double[n];
	int i;
#	pragma omp parallel for num_threads(OMP_THREAD_NUM) \
	default(none) private(i) shared(Bx, ones, pc, Y)
	for(i = 0; i < n; i++){
		Bx[i] = 0; Y[i] = pc[0];
		ones[i] = 1;
	}
	//Bx = c1 * X
	cblas_daxpy(n, c1, X, 1, Bx, 1);
	//Bx = Bx + c0
	cblas_daxpy(n, c0, ones, 1, Bx, 1);
	for(i = 0; i < deg; i++){
		//Y = Bx .* Y
		vdMul(n, Bx, Y, Y);
		//Y = Y + pc[i+1]
		cblas_daxpy(n, pc[i+1], ones, 1, Y, 1);
	}
	delete []Bx;
	delete []ones;
}

/* ployChebyValue - compute chebyshev polynomial p(X), p() is applied to each element of X
 *   Inputs :
 *         n    - the size of X
 *         X    - array, size at least n
 *         a, b - interval [a, b]
 *         deg  - the degree of the chebyshev polynomial
 *   Outputs:
 *         Y    - preallocated array, size at least n
 */
void polyChebyValue(const MKL_INT n, const double* X, double a, double b, const MKL_INT deg, double* Y){
	double c0 = (a + b) / (a - b), c1 = 2 / (b - a);
	double *Bx;
	double *ones;
	double *T1, *T2;
	int i;
#	pragma omp parallel for num_threads(OMP_THREAD_NUM) \
	default(none) private(i) shared(Y)
	for(i = 0; i < n; i++){
		Y[i] = 1;
	}
	if(deg == 0){
		return;
	}
	if(deg == 1){
		cblas_dscal(n, c0, Y, 1);
		cblas_daxpy(n, c1, X, 1, Y, 1);
		return;
	}
	// deg >= 2
	Bx = new double[n];
	T1 = new double[n];
	T2 = new double[n];
#	pragma omp parallel for num_threads(OMP_THREAD_NUM) \
	default(none) private(i) shared(T1, Bx, c0)
	for(i = 0; i < n; i++){
		T1[i] = 1;
		Bx[i] = c0;
	}
	//Bx := Bx + c1 * X
	cblas_daxpy(n, c1, X, 1, Bx, 1);
	//T2 := Bx
	cblas_dcopy(n, Bx, 1, T2, 1);
	//Bx := 2*Bx
	cblas_dscal(n, 2, Bx, 1);
	for(i = 2; i <= deg; i++){
		if(i % 2 == 0){
			//Y := 2Bx * T2
			vdMul(n, Bx, T2, Y);
			//Y := Y - T1
			cblas_daxpy(n, -1, T1, 1, Y, 1);
			if(i < deg)	cblas_dcopy(n, Y, 1, T1, 1); // T1 := Y
		}
		else{
			//Y := 2Bx * T1
			vdMul(n, Bx, T1, Y);
			//Y := Y - T2
			cblas_daxpy(n, -1, T2, 1, Y, 1);
			if(i < deg) cblas_dcopy(n, Y, 1, T2, 1); // T2 := Y
		}
	}
	delete []Bx;
	delete []T1;
	delete []T2;
}

/* update_degree - update polynomial degree
 *   Input :
 *        lk, lkw - best eigenvalue, lk is the k-th largest eigenvalue
 *        a, b    - interval [a, b]
 *        pcs     - polynomial coefficients, p[i] is an array, when using chebyshev polynomials, set this parameter to NULL
 *        maxdeg  - max degree of the polynomial
 *   Output :
 *        deg     - updated degree
 */
void update_degree(double lk, double lkw, double a, double b, const double** pcs, const MKL_INT maxdeg, int* deg){
	int degm;
	double evs[] = {lk, lkw}, ev_res[2];
	for(degm = 2;degm <= maxdeg; degm++){
		if(pcs == NULL){ //using chebyshev
			polyChebyValue(2, evs, a, b, degm, ev_res);
		}
		else{ //using general polynomials
			polyValue(2, evs, a, b, pcs[deg], degm, ev_res);
		}
		//check breaking rules
		if(ev_res[2] / ev_res[1] < 0.9) break;
	}
	//output
	*deg = degm > maxdeg ? maxdeg : degm;
}
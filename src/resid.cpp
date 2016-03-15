#include "arrabit_inc.h"

/* residual - given Z, compute relative residuals, i.e. AZ-mZ / m 
 *   Inputs :
 *         n, k                    - A is n-by-n, Z is n-by-k
 *         val, indx, pntrB, pntrE - represents the sparse matrix A
 *         Z                       - array, size at least n-by-k
 *         d                       - array, size at least k, stores all eigenvalues
 *   Outputs :
 *         res                     - pre-allocated array, size at least k
 */
void residual(const MKL_INT n, const MKL_INT k, const double* var, const MKL_INT* indx, const MKL_INT* pntrB, const MKL_INT* pntrE, const double *Z, const double *d, double* res){
	char matdescra[] = {'G', 'L', 'N', 'F'};
	char trans[] = {'N'};
	double alpha = 1.0, beta = 0.0;
	int i;
	double AZ = new double[n * k];
	//AZ = A * Z
	mkl_dcscmm(trans, &n, &k, &n, &alpha, matdescra, val, indx, pntrB, pntrE, Z, &n, &beta, AZ, &n);
	//AZ - d[i] * Z
#	pragma omp parallel for num_threads(OMP_THREAD_NUM) \
	default(none) private(i) shared(AZ, Z)
	for(i = 0; i < k; i++){
		cblas_daxpy(n, -d[i], Z, 1, AZ, 1);
	}
	//relative residual
#	pragma omp parallel for num_threads(OMP_THREAD_NUM) \
	default(none) private(i) shared(AZ)
	for(i = 0; i < n; i++){
		res[i] = cblas_dnrm2(n, AZ + i * n, 1);
		// res[i] := res[i] / max(1, abs(d[i]))
		res[i] /= (fabs(d[i]) > 1 ? fabs(d[i]) : 1);
		}
	}
}
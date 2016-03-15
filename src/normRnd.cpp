#include "mkl.h"
#include <time.h>

/* normRnd - generate random Gaussian matrix, feedback shift generator, box-muller method
 *   Inputs:
 *         m, n  - size of the matrix
 *         X     - pre-allocated array, size at least m-by-n
 *         mu    - expectation
 *         sigma - std. error
 *   Outputs:
 *         X     - random matrix
 */
void normRnd(const MKL_INT m, const MKL_INT n, double* X, double mu, double sigma){
	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_R250, (MKL_UINT)time(NULL));
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, m * n, X, mu, sigma);
	vslDeleteStream(&stream);
}
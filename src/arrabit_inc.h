#include "mkl.h"
#include "spmat.h"
#include <omp.h>
#include <math.h>
#include <time.h>

/* Developer tags : file|code|debug
 *   file: source filename
 *   code: Y/N
 *   debug: Y/N
 *---------------------
 * Eg: function 'normRnd' is implemented in normRnd.cpp, the coding work is done, but the debugging work hasn't started
 * should go like this:
 *   normRnd.cpp|Y|N
 */

extern int OMP_THREAD_NUM;
//extern functions
//eigs_sps.cpp|N|N
extern "C" int eigs_sps(MKL_INT, double*, MKL_INT*, MKL_INT*, MKL_INT*, char*, MKL_INT, double*);

//arrabit.cpp|N|N
int arrabit(const MKL_INT, const MKL_INT, const double*, const MKL_INT*, const MKL_INT*, const MKL_INT*, const char, Opts*, double*, double*);
void call_innersolver();

//normRnd.cpp|Y|Y
void normRnd(const MKL_INT, const MKL_INT, double*, double, double);

//arrabit.cpp|Y|Y
void mat_dcopy(const MKL_INT, const MKL_INT, const double *dest, double *target);

//norm.cpp|Y|N
double norm1(double*, const MKL_INT, const MKL_INT);
void normCol2(double*, const MKL_INT, const MKL_INT);

//poly.cpp|Y|N
void polyAX(const MKL_INT, const MKL_INT, const double*, const MKL_INT*, const MKL_INT*, const MKL_INT*, const double*, double, double, const double*, const MKL_INT, double*);
void polyChebyAX(const MKL_INT, const MKL_INT, const double*, const MKL_INT*, const MKL_INT*, const MKL_INT*, const double*, double, double, const MKL_INT, double*);
void polyValue(const MKL_INT, const double*, double, double, const double*, const MKL_INT, double*);
void polyChebyValue(const MKL_INT, const double*, double, double, const MKL_INT, double*);
void update_degree(double, double, double, double, const double**, const MKL_INT, int*);

//aug.cpp|Y|N
void doAug(const MKL_INT, const MKL_INT, const double*, const MKL_INT*, const MKL_INT*, const MKL_INT*, cosnt double*, const int, double, double*);
void proj2Qker(const MKL_INT, const MKL_INT, const MKL_INT, const MKL_INT, double*, const double*, const double*);

//RR.cpp|Y|N
void RR(const MKL_INT, const MKL_INT, const double*, const MKL_INT*, const MKL_INT*, const MKL_INT*, double*, double*, double*);

//eigqr.cpp|Y|Y
void eigqr(const MKL_INT, double*, double*, double*);

//resid.cpp|Y|N
void residual(const MKL_INT, const MKL_INT, const double*, const MKL_INT*, const MKL_INT*, const MKL_INT*, const double *, double*);

//est.cpp|Y|N
void est_key_eigenvalues(const MKL_INT, const MKL_INT, double*, MKL_INT*, MKL_INT*, MKL_INT*, double*, int, int,  double*, double*, double*);

typedef struct {
	int maxit0;
	int maxit1;
	int maxit2;
	int cheby;
	int deg;
	int maxdeg;
	int doARR;
	int fixdeg;
	int nblks;
	int maxblks;
	double guard;
	char inner_solver;
	double tol;
	double tol_t;
	double tiny;
	double tiny1;
	double tiny2;
} Opts;
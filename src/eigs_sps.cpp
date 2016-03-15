#include<mkl.h>
#include<stdio.h>

//extern fortran subroutines
//call dsaupd  
//     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
//      IPNTR, WORKD, WORKL, LWORKL, INFO )
//
extern "C" void dsaupd_(int*, char*, int*, char*, int*, double*, double*, int*, double*, int*, int*, int*, double*, double*, int*, int*, int, int);

//extern fortran subroutines
//  call dseupd  
//    ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
//      RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
extern "C" void dseupd_(int*, char*, int*, double*, double*, int*, double*, char*, int*, char*, int*, double*, double*, int*, double*, int*, int*, int*, double*, double*, int*, int*, int, int, int);

extern "C" int eigs_sps(MKL_INT, double*, MKL_INT*, MKL_INT*, MKL_INT*, char*, MKL_INT, double*);
/* eigs_sps - computs a few eigenvalues of a sparse symmetric matrix
 *   Inputs :
 *         n                       - size
 *         val, indx, pntrB, pntrE - represents the sparse matrix A
 *         which                   - character array, can be 'SA' 'LA' 'BE' 'SM' 'LM'
 *         nev                     - number of eigenvalues to be computed
 *         d                       - pre-allocated array, serves as the output eigenvalues
 *   Return Value :
 *         (int) 0                 - process exits normally
 *         (int) != 0              - process exits with errors
 */
int eigs_sps(MKL_INT n, double* val, MKL_INT* indx, MKL_INT* pntrB, MKL_INT* pntrE, char* which, MKL_INT nev, double *d){
	int ido = 0, ncv = int(nev*1.1+1), ldv = n, lworkl = ncv*(ncv+8), info = 0, revc = 0, ierr = 0;
	int iparam[11], ipntr[11], *select;
	unsigned long pntr0, pntr1;
	double tol = 0, alpha = 1.0, beta = 0.0, sigma = 0.0;
	double *resid, *workd, *workl, *v;
	char bmat[] = {'I'}, transa[] = {'N'}, matdescra[] = {'G', 'L', 'N', 'F'}, howmny[] = {'A'};
	resid = new double[n];
	workd = new double[3*n];
	workl = new double[lworkl];
	v = new double[ncv*ldv];
	select = new int[ncv];
	iparam[0] = 1;
	iparam[2] = 300;
	iparam[6] = 1;
	while(true){
		dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info, 1, 2);
		if(ido == 1 || ido == -1){
			pntr0 = (unsigned long)ipntr[0]; pntr1 = (unsigned long)ipntr[1];
			mkl_dcscmv(transa, &n, &n, &alpha, matdescra, val, indx, pntrB, pntrE, (double*)pntr0, &beta, (double*)pntr1);
		}
	}
	dseupd_(&revc, howmny, select, d, v, &ldv, &sigma, bmat, &n, which, &ncv, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &ierr, 1, 1, 2);
	if(ierr != 0){
		printf("Error with dseupd, info = %d\n", ierr);
	}
	delete []resid;
	delete []workd;
	delete []workl;
	delete []select;
	delete []v;
	return info;
}

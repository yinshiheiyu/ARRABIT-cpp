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

extern "C" int eigs_ges(MKL_INT, double*, MKL_INT, char*, MKL_INT, double*);
/* eigs_sps - computs a few eigenvalues of a general symmetric matrix
 *   Inputs :
 *         n     - size
 *         A     - array, represents matrix A, ColMajor
 *         lda   - leading dimention of A
 *         which - character array, can be 'SA' 'LA' 'BE' 'SM' 'LM'
 *         nev   - number of eigenvalues to be computed
 *         d     - pre-allocated array, serves as the output eigenvalues
 *   Return Value :
 *         (int) 0                 - process exits normally
 *         (int) != 0              - process exits with errors
 */
int eigs_ges(MKL_INT n, double* A, MKL_INT lda, char* which, MKL_INT nev, double *d){
	int ido = 0, ncv = int(nev*1.1+1), ldv = n, lworkl = ncv*(ncv+8), info = 0, revc = 0, ierr = 0;
	int iparam[11], ipntr[11], *select;
	double tol = 0, alpha = 1.0, beta = 0.0, sigma = 0.0;
	double *resid, *workd, *workl, *v;
	char bmat[] = {'I'}, howmny[] = {'A'};
	resid = new double[n];
	workd = new double[3*n];
	workl = new double[lworkl];
	v = new double[ncv*ldv];
	select = new int[ncv];
	iparam[0] = 1;
	iparam[2] = 300;
	iparam[6] = 1;
	//ensure ncv .le. n
	ncv = ncv < n ? ncv : n;
	while(true){
		dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info, 1, 2);
		if(ido == 1 || ido == -1){
			printf("%d, %d\n", ipntr[0], ipntr[1]);
			cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, A, lda, workd + ipntr[0] - 1, 1, beta, workd + ipntr[1] - 1, 1);
			printf("cblas_dgemv subroutine done.\n");
		}
		if(ido == 99) break;
	}
	printf("dsaupd_ subroutine done.\n");
	printf("%d, %d, %d\n", nev, ncv, n);
	dseupd_(&revc, howmny, select, d, v, &ldv, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &ierr, 1, 1, 2);
	printf("dseupd_ subroutine done.\n");
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

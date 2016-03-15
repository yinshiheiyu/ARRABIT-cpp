typedef struct{
	char *matrix_filename;
	char *result_filename;
	int maxit0, maxit1, maxit2;
	int nblks; //augmentation blocks
	int maxblks;
	int nev; //number of eigenvalues, i.e. k
	double guard; //rate of guard vectors, q = 0.1k
	int doARR; //whether to do augmentation
	int defl; //whether to do deflation
	int trim; //whether to do trim
	int deg0; //initial degree
	int maxdeg; //max degree
	double tol;
	double tol0;
	int tolint;
	double *l1;
	double *ln;
	double *lkw;
	char innersolver; //indicated by char. 'M' for MPM and 'G' for GN
	//other fields
} OPTS;
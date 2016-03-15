#include<stdio.h>
/* dsort - a wrapper function for dsortr in ARPACK
 * external funtion  :  dsortr( WHICH, APPLY, N, X1, X2 )
 *   Inputs :
 *         which       - characters, can be 'SA' 'LA' 'SM' 'LM'
 *         permutation - specifies whether permutation array is required (0/1)
 *         n           - size of the array
 *         x1          - the array to be sorted
 *         x2          - the permutation array, when permutation .eq. 0, it is not referenced
 */

extern "C" void dsortr_(char*, int*, int*, double*, double*, int);

extern "C" void dsort(char*, int, int, double*, double*);
/*
int main(){
	double x1[] = {3.2, 2.2, 5.6, 7, 8.8, 1.0, 0.4};
	double x2[7];
	char which[] = {'S', 'A'};
	for(int i = 0; i < 7; i++) printf("%f\t", x1[i]);
	printf("\n");
	dsort(which, 1, 7, x1, x2);
	for(int i = 0; i < 7; i++) printf("%f\t", x1[i]);
	printf("\n");
	for(int i = 0; i < 7; i++) printf("%f\t", x2[i]);
	printf("\n");
	return 0;
}*/

void dsort(char* which, int permutation, int n, double* x1, double* x2){
	for(int i = 0; i < n; i++) x2[i] = i;
	dsortr_(which, &permutation, &n, x1, x2, 2);
}

#define min(x,y) (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) < (y)) ? (y) : (x))
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "opts.h"


void read_opts(const char *opts_filename, OPTS *opts);
OPTS* createEmptyOpts(void);

int main(){
	char *opts_filename = "parameters.in";
	OPTS *opts;
	opts = createEmptyOpts();
	read_opts(opts_filename, opts);
	//initialize
	opts->tol0 = max(opts->tol, 1e-6);
	opts->maxblks = min(4, 2 * opts->nblks);

	printf("OPTS File has been set.\n");
	printf("%c\n", opts->innersolver);
	printf("%s\n", opts->matrix_filename);
	printf("%lf\n", opts->tol);
	printf("%lf\n", opts->tol0);
	if (opts->ln){
		printf("%lf\n", *opts->ln);
	}
	else{
		printf("ln should be estimate!\n");
	}
	printf("%d\n", opts->maxit0);
	getchar();
	return 0;
}

OPTS* createEmptyOpts(void){
	OPTS *opts = (OPTS *)malloc(sizeof(OPTS));
	if (opts != NULL)
	{
		//initialize
		opts->matrix_filename = NULL;
		opts->result_filename = NULL;
		opts->maxit0 = 30;
		opts->maxit1 = 10;
		opts->maxit2 = 5;
		opts->nblks = 2;
		opts->maxblks = min(4,2*opts->nblks);
		opts->guard = 0.1;
		opts->nev = 10;
		opts->doARR = 1;
		opts->defl = 1;
		opts->trim = 1;
		opts->tol = 1e-6;
		opts->deg0 = 3;
		opts->maxdeg = 15;
		opts->tol = 1e-6;
		opts->tol0 = max(opts->tol, 1e-6);
		opts->tolint = 2;
		opts->l1 = NULL;
		opts->ln = NULL;
		opts->lkw = NULL;
		opts->innersolver = 'M';
		return opts;
	}
	printf("Out of space!!\n");
	return NULL;
}

void read_opts(const char *opts_filename, OPTS *opts){
	FILE *fp;
	char line[256];
	char line_var[256];
	char line_value[256];
	int nev = 0;

	if ((fp = fopen(opts_filename, "r")) != NULL){
		while (fgets(line, 256, fp) != NULL){
			sscanf(line, "%s %s", line_var, line_value);
			
			if (!strcmp(line_var, "matrix_filename")){
				opts->matrix_filename = (char *)malloc(sizeof(char));
				if (opts->matrix_filename){
					strcpy(opts->matrix_filename, line_value);
					printf("%s\n", opts->matrix_filename);
				}
				else{
					printf("Out of space!!\n");
					break;
				}
			}
			else if (!strcmp(line_var, "result_filename")){
				opts->result_filename = (char *)malloc(sizeof(char));
				if (opts->result_filename){
					strcpy(opts->result_filename, line_value);
					printf("%s\n", opts->result_filename);
				}
				else{
					printf("Out of space!!\n");
					break;
				}
			}
			else if (!strcmp(line_var, "maxit0")){
				sscanf(line_value, "%d", &opts->maxit0);
				printf("%d\n", opts->maxit0);
			}
			else if (!strcmp(line_var, "maxit1")){
				sscanf(line_value, "%d", &opts->maxit1);
				printf("%d\n", opts->maxit1);
			}
			else if (!strcmp(line_var, "maxit2")){
				sscanf(line_value, "%d", &opts->maxit2);
				printf("%d\n", opts->maxit2);
			}
			else if (!strcmp(line_var, "nblks")){
				sscanf(line_value, "%d", &opts->nblks);
				printf("%d\n", opts->nblks);
			}
			else if (!strcmp(line_var, "maxblks")){
				sscanf(line_value, "%d", &opts->maxblks);
				printf("%d\n", opts->maxblks);
			}
			else if (!strcmp(line_var, "nev")){
				sscanf(line_value, "%d", &opts->nev);
				printf("%d\n", opts->nev);
			}
			else if (!strcmp(line_var, "guard")){
				sscanf(line_value, "%lf", &opts->guard);
				printf("%lf\n", opts->guard);
			}
			else if (!strcmp(line_var, "doARR")){
				sscanf(line_value, "%d", &opts->doARR);
				printf("%d\n", opts->doARR);
			}
			else if (!strcmp(line_var, "defl")){
				sscanf(line_value, "%d", &opts->defl);
				printf("%d\n", opts->defl);
			}
			else if (!strcmp(line_var, "trim")){
				sscanf(line_value, "%d", &opts->trim);
				printf("%d\n", opts->trim);
			}
			else if (!strcmp(line_var, "deg0")){
				sscanf(line_value, "%d", &opts->deg0);
				printf("%d\n", opts->deg0);
			}
			else if (!strcmp(line_var, "maxdeg")){
				sscanf(line_value, "%d", &opts->maxdeg);
				printf("%d\n", opts->maxdeg);
			}
			else if (!strcmp(line_var, "tol")){
				sscanf(line_value, "%lf", &opts->tol);
				printf("%lf\n", opts->tol);
			}
			else if (!strcmp(line_var, "tol0")){
				sscanf(line_value, "%lf", &opts->tol0);
				printf("%lf\n", opts->tol0);
			}
			else if (!strcmp(line_var, "tolint")){
				sscanf(line_value, "%d", &opts->tolint);
				printf("%d\n", opts->tolint);
			}
			else if (!strcmp(line_var, "l1")){
				double l1;
				opts->l1 = (double *)malloc(sizeof(double));
				if (opts->l1){
					sscanf(line_value, "%lf", &l1);
					*opts->l1 = l1;
					printf("%lf\n", opts->l1);
				}
				else{
					printf("Out of space!!\n");
					break;
				}
			}
			else if (!strcmp(line_var, "ln")){
				double ln;
				opts->ln = (double *)malloc(sizeof(double));
				if (opts->ln){
					sscanf(line_value, "%lf", &ln);
					*opts->ln = ln;
					printf("%lf\n", opts->ln);
				}
				else{
					printf("Out of space!!\n");
					break;
				}
			}
			else if (!strcmp(line_var, "lkw")){
				double lkw;
				opts->lkw = (double *)malloc(sizeof(double));
				if (opts->lkw){
					sscanf(line_value, "%lf", &lkw);
					*opts->lkw = lkw;
					printf("%lf\n", opts->lkw);
				}
				else{
					printf("Out of space!!\n");
					break;
				}
			}
			else if (!strcmp(line_var, "innersolver")){
				sscanf(line_value, "%c", &opts->innersolver);
				printf("%c\n", opts->innersolver);
			}
		}
		fclose(fp);
	}
	else{
		printf("ERROR: No such File named \"%s\"", opts_filename);
	}
}
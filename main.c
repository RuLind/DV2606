#include <stdio.h>
#include <pthread.h>
#include <string.h>
#include <stdlib.h>
#define MAX_SIZE 4096

typedef double matrix[MAX_SIZE][MAX_SIZE];

int	N;		/* matrix size		*/
int	maxnum;		/* max number of element*/
char	*Init;		/* matrix init type	*/
int	PRINT;		/* print switch		*/
matrix	A;		/* matrix A		*/
double	b[MAX_SIZE];	/* vector b             */
double	y[MAX_SIZE];	/* vector y             */

/* forward declarations */
void work(void);
void Init_Matrix(void);
void Print_Matrix(void);
void Init_Default(void);
int Read_Options(int, char **);
struct threadArgs {
    unsigned int comp_row;
    unsigned int start_col;
    unsigned int nr_of_col;
};
int
main(int argc, char **argv)
{
    int i, timestart, timeend, iter;

    Init_Default();		/* Init default values	*/
    Read_Options(argc,argv);	/* Read arguments	*/
    Init_Matrix();		/* Init the matrix	*/
    work();
    if (PRINT == 1)
	   Print_Matrix();
}
void* child_labour(void* params){
    int i,j, start_col, nr_of_col,comp_row;
    struct threadArgs *args = (struct threadArgs*) params;
    start_col = args->start_col;
    nr_of_col = args->nr_of_col;
    comp_row = args->comp_row;
    for (i = start_col; i < start_col+nr_of_col; i++) {
	        for (j = comp_row+1; j < N; j++)
		        A[i][j] = A[i][j] - A[i][comp_row]*A[comp_row][j]; /* Elimination step */
	        b[i] = b[i] - A[i][comp_row]*y[comp_row];
	        A[i][comp_row] = 0.0;
	    }
}
void
work(void)
{
    int i, j, k, rest, numThreads, row_per_thread;
    pthread_t* children;
    struct threadArgs* args;
    children = malloc(numThreads * sizeof(pthread_t));
    numThreads = 16;//matching number of threads on cpu
    // children = malloc(numThreads * sizeof(pthread_t));
    /* Gaussian elimination algorithm, Algo 8.4 from Grama */
    for (k = 0; k < N; k++) { /* Outer loop */
        if(k+16>N){
            numThreads=numThreads-1;
            //printf("we have %d rows and %d threads\n",N-k,numThreads);
        }
	    for (j = k+1; j < N; j++)
	       A[k][j] = A[k][j] / A[k][k]; /* Division step */
	    y[k] = b[k] / A[k][k];
	    A[k][k] = 1.0;
        row_per_thread = (N-k)/numThreads; //heltals division
        rest = (N-k)%numThreads;
        int test = row_per_thread*numThreads +rest;
        for(i = 0; i<numThreads; i++){ //split sub step into numThreads
            args = malloc(sizeof(struct threadArgs));
            args->comp_row = k;
            args->start_col = i*row_per_thread+k+1;
            args->nr_of_col = row_per_thread;
            if(i==0){
                args->nr_of_col = row_per_thread + rest;
            }else{
                args->start_col = args->start_col+rest;
            }
            //printf("thread %d is between %d to %d, and has %d rows to go through\n",i,args->start_col,args->start_col+args->nr_of_col,args->nr_of_col);
            pthread_create(&(children[i]), // our handle for the child
                NULL, // attributes of the child
                child_labour, // the function it should run
                (void*)args);
        }
        }
        for (unsigned int id = 0; id < numThreads; id++) {
		    pthread_join(children[id], NULL );
        }
        free(children);
        free(args);
    }


void
Init_Matrix()
{
    int i, j;

    printf("\nsize      = %dx%d ", N, N);
    printf("\nmaxnum    = %d \n", maxnum);
    printf("Init	  = %s \n", Init);
    printf("Initializing matrix...");

    if (strcmp(Init,"rand") == 0) {
        for (i = 0; i < N; i++){
            for (j = 0; j < N; j++) {
                if (i == j) /* diagonal dominance */
                    A[i][j] = (double)(rand() % maxnum) + 5.0;
                else
                    A[i][j] = (double)(rand() % maxnum) + 1.0;
            }
        }
    }
    if (strcmp(Init,"fast") == 0) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (i == j) /* diagonal dominance */
                    A[i][j] = 5.0;
                else
                    A[i][j] = 2.0;
            }
        }
    }

    /* Initialize vectors b and y */
    for (i = 0; i < N; i++) {
        b[i] = 2.0;
        y[i] = 1.0;
    }

    printf("done \n\n");
    if (PRINT == 1)
        Print_Matrix();
}

void
Print_Matrix()
{
    int i, j;

    printf("Matrix A:\n");
    for (i = 0; i < N; i++) {
        printf("[");
        for (j = 0; j < N; j++)
            printf(" %5.2f,", A[i][j]);
        printf("]\n");
    }
    printf("Vector b:\n[");
    for (j = 0; j < N; j++)
        printf(" %5.2f,", b[j]);
    printf("]\n");
    printf("Vector y:\n[");
    for (j = 0; j < N; j++)
        printf(" %5.2f,", y[j]);
    printf("]\n");
    printf("\n\n");
}

void
Init_Default()
{
    N = 2048;
    Init = "rand";
    maxnum = 15.0;
    PRINT = 0;
}

int
Read_Options(int argc, char **argv)
{
    char    *prog;

    prog = *argv;
    while (++argv, --argc > 0)
        if (**argv == '-')
            switch ( *++*argv ) {
                case 'n':
                    --argc;
                    N = atoi(*++argv);
                    break;
                case 'h':
                    printf("\nHELP: try sor -u \n\n");
                    exit(0);
                    break;
                case 'u':
                    printf("\nUsage: gaussian [-n problemsize]\n");
                    printf("           [-D] show default values \n");
                    printf("           [-h] help \n");
                    printf("           [-I init_type] fast/rand \n");
                    printf("           [-m maxnum] max random no \n");
                    printf("           [-P print_switch] 0/1 \n");
                    exit(0);
                    break;
                case 'D':
                    printf("\nDefault:  n         = %d ", N);
                    printf("\n          Init      = rand" );
                    printf("\n          maxnum    = 5 ");
                    printf("\n          P         = 0 \n\n");
                    exit(0);
                    break;
                case 'I':
                    --argc;
                    Init = *++argv;
                    break;
                case 'm':
                    --argc;
                    maxnum = atoi(*++argv);
                    break;
                case 'P':
                    --argc;
                    PRINT = atoi(*++argv);
                    break;
                default:
                    printf("%s: ignored option: -%s\n", prog, *argv);
                    printf("HELP: try %s -u \n\n", prog);
                    break;
            }
}

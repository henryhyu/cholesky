#include "csparse.h"
#include "mmio.h"
#include <stdlib.h>
# include <stdio.h>
#include <time.h>
int main (int argc, char *argv[]) {
    int ret_code;
	FILE *f;
	int row, col, nz;  
    int i, *I, *J;
    double *val, *b;
   // cs *tmp, *A;

    // Read in input mtx
    if ((f = fopen(argv[1], "r")) == NULL)
    	exit(1);

    if ((ret_code = mm_read_mtx_crd_size(f, &row, &col, &nz)) != 0)
        exit(1);

    I = (int *)malloc(nz * sizeof(int));
    J = (int *)malloc(nz * sizeof(int));
    val = (double *)malloc(nz * sizeof(double));

    for (i = 0; i < nz; i++) {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;
        J[i]--;
    }
    fclose(f);

    // set up matrix
    cs *tmp = cs_spalloc(row, col, nz, 1, 1);
    for (i = 0; i < nz; i++) {
      cs_entry(tmp, I[i], J[i], val[i]);
    }
    cs *A = cs_triplet(tmp);
    int n = A->n ;
    
    // set up b
    double *realx = cs_calloc(n, sizeof(double));
    realx[0]=-1;
    realx[row-1]=1;
    //b = (double *)malloc(row * sizeof(double));
    b = cs_calloc(n, sizeof(double));
    //for (i = 0; i < row; i++) b[i] = 0;

    cs_gaxpy(A,realx,b);
    
    
    
    // prints A
    // cs_print(A, 0);

    // cholesky solve
    // printf("ok: %d\n", cs_cholsol(A, b, row));
    // for (i = 0; i < row; i++) printf("%lg ", b[i]);
    // printf("\n");
    /*
    EXPECTING: 
    x1 = 1.5
    x2 = 2
    x3 = 1.5
    */
    // return 0;
    // broken down cholesky
    double *x ;
    css *S ;
    csn *N ;
    int ok, order;
    order = row;
    if (!A || !b) return (1) ;      /* check inputs */

    clock_t symb_time, num_time, solve_time, total_time;
    clock_t start = clock();
    clock_t start2 = start;
    S = cs_schol (A, order) ;       /* ordering and symbolic analysis */
    symb_time = clock()-start;
    start = clock();
    N = cs_chol (A, S) ;        /* numeric Cholesky factorization */
    num_time = clock()-start;
    x = cs_malloc (n, sizeof (double));
    ok = (S && N && x) ;
    if (ok)
    {
      start = clock();
      cs_ipvec (n, S->Pinv, b, x) ;   /* x = P*b */
      cs_lsolve (N->L, x) ;       /* x = L\x */
      cs_ltsolve (N->L, x) ;      /* x = L'\x */
      cs_pvec (n, S->Pinv, x, b) ;    /* b = P'*x */
      solve_time = clock()-start;
    }
    total_time = clock()-start2;
    printf("%d\n", ok);
    printf("Time taken on ordering and symbolic analysis: %f s\n", symb_time*1.0/CLOCKS_PER_SEC);
    printf("Time taken on numeric factorization: %f s\n", num_time*1.0/CLOCKS_PER_SEC);
    printf("Time taken on triangular solve: %f s\n", solve_time*1.0/CLOCKS_PER_SEC);
    printf("Total time taken: %f s\n", total_time*1.0/CLOCKS_PER_SEC);

    double error;
    for (i = 0; i < row; i++) {
      error+=pow(b[i]-realx[i],2);
    }
    printf("Error in lhs: %f\n", sqrt(error));
    
    cs_free (realx) ;
    cs_free (x) ;
    cs_sfree (S) ;
    cs_nfree (N) ;
    
    
    return 0;
}

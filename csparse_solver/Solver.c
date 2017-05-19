#include "csparse.h"
#include "mmio.h"
#include <stdlib.h>
# include <stdio.h>

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

	// set up b
    b = (double *)malloc(row * sizeof(double));
    for (i = 0; i < row; i++) b[i] = 1;
    
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
    int n, ok, order;
    order = row;
    if (!A || !b) return (1) ;      /* check inputs */
    n = A->n ;
    S = cs_schol (A, order) ;       /* ordering and symbolic analysis */
    N = cs_chol (A, S) ;        /* numeric Cholesky factorization */
    x = cs_malloc (n, sizeof (double));
    ok = (S && N && x) ;
    if (ok)
    {
    cs_ipvec (n, S->Pinv, b, x) ;   /* x = P*b */
    cs_lsolve (N->L, x) ;       /* x = L\x */
    cs_ltsolve (N->L, x) ;      /* x = L'\x */
    cs_pvec (n, S->Pinv, x, b) ;    /* b = P'*x */
    }
    cs_free (x) ;
    cs_sfree (S) ;
    cs_nfree (N) ;
    printf("%d\n", ok);
    return 0;
}

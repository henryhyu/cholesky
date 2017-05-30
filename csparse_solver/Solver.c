#include "csparse.h"
#include "mmio.h"
#include <stdlib.h>
# include <stdio.h>
#include <time.h>
int main (int argc, char *argv[]) {
    int ret_code;
    FILE *f;
    int row, col, nz;  
    int i, *I, *J, n;
    double *val, *b;
    int order=0;

    if(argc > 2) {
      order = atoi(argv[2]);
    }

    // Read in input mtx
    if ((f = fopen(argv[1], "r")) == NULL)
        exit(1);

    const char *dot = strrchr(argv[1], '.');
    if (!strcmp(dot,".mtx")) {
      if ((ret_code = mm_read_mtx_crd_size(f, &row, &col, &nz)) != 0)
        exit(1);
    }
    else if (!strcmp(dot,".bin")) {
      int tempn;
      fread(&tempn, sizeof(int), 1, f);
      fread(&nz, sizeof(int), 1, f);
      row = tempn;
      col = tempn;
    }

    I = (int *)malloc(nz * sizeof(int));
    J = (int *)malloc(nz * sizeof(int));
    val = (double *)malloc(nz * sizeof(double));
    
    if (!strcmp(dot,".mtx")) {
      for (i = 0; i < nz; i++) {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;
        J[i]--;
      }
    }
    else if(!strcmp(dot,".bin")) {
      for (i = 0; i < nz; i++) {
        fread(&I[i], sizeof(int), 1, f);
        fread(&J[i], sizeof(int), 1, f);
        fread(&val[i], sizeof(double), 1, f);
        I[i]--;
        J[i]--;
      }
    }

    fclose(f);
    
    // set up triplet matrix
    cs *tmp = cs_spalloc(row, col, nz, 1, 1);
    for (i = 0; i < nz; i++) {
      cs_entry(tmp, I[i], J[i], val[i]);
    }
    cs *A = cs_triplet(tmp);
    A->nz = tmp->nz;
    n = A->n ;
    
    // set up b
    double *realx = cs_calloc(n, sizeof(double));
    realx[0]=-1;
    realx[n-1]=1;
    b = cs_calloc(n, sizeof(double));
    cs_gaxpy(A,realx,b);
    
    // broken down cholesky
    double *x ;
    css *S ;
    csn *N ;
    int ok;

    if (!A || !b) return (1);      /* check inputs */
    
    clock_t symb_time, num_time, solve_time, total_time;
    clock_t start = clock();
    clock_t start2 = start;
    S = cs_schol (A, order) ;       /* ordering and symbolic analysis */
    symb_time = clock()-start;
    start = clock();
    N = cs_chol (A, S);        /* numeric Cholesky factorization */
    num_time = clock()-start;
    x = cs_malloc (n, sizeof (double));
    ok = (S && N && x);
//            printf("hello %p,%p,%p\n", S, N, x);
    fflush(stdout);
    if (ok) {
      start = clock();
      cs_ipvec (n, S->Pinv, b, x) ;   /* x = P*b */
      cs_lsolve (N->L, x) ;       /* x = L\x */
      cs_ltsolve (N->L, x) ;      /* x = L'\x */
      cs_pvec (n, S->Pinv, x, b) ;    /* b = P'*x */
      solve_time = clock()-start;
    }
    total_time = clock()-start2;
//   printf("%d\n", ok);
    printf("%s,", argv[1]);
   printf("%f,", symb_time*1.0/CLOCKS_PER_SEC); // Time taken on ordering and symbolic analysis
   printf("%f,", num_time*1.0/CLOCKS_PER_SEC); // Time taken on numeric factorization
   printf("%f,", solve_time*1.0/CLOCKS_PER_SEC); // Time taken on triangular solve
   printf("%f,", total_time*1.0/CLOCKS_PER_SEC); // Total time taken: 
   // double error;
   // for (i = 0; i < n; i++) {
   //   error+=pow(b[i]-realx[i],2);
   // }
   // printf("Error in lhs: %f\n", sqrt(error));
    
    printf("%d,", A->nzmax); //nnz(A)
    printf("%d\n", N->L->nzmax); // nnz(R)
    //printf("nnz(R) / nnz(A): %f\n", 1.0 * N->L->nzmax / A->nzmax);

    cs_free (realx) ;
    cs_free (x) ;
    cs_free (b) ;
    cs_sfree (S) ;
    cs_nfree (N) ;
    cs_spfree (A);
    // can be freed earlier if need memory
    cs_spfree (tmp);
    free(I);
    free(J);
    free(val);

    return 0;
}
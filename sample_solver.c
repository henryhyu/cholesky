#include <stdio.h>
#include <stdlib.h>
#include "csparse.h"
#include "mmio.h"

int main (int argc, char *argv[]) {
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nnz;   
    int i, *I, *J;
    double *val;
    
    if (argc < 2) {
      fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
      exit(1);
    }
    else { 
      if ((f = fopen(argv[1], "r")) == NULL) 
	exit(1);
    }
    
    if (mm_read_banner(f, &matcode) != 0) {
      printf("Could not process Matrix Market banner.\n");
      exit(1);
    }
    
    
    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
	mm_is_sparse(matcode) ) {
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(1);
    }
    
    /* find out size of sparse matrix .... */
    
    /* reseve memory for matrices */
    
    I = (int *) malloc(nnz * sizeof(int));
    J = (int *) malloc(nnz * sizeof(int));
    val = (double *) malloc(nnz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i < nnz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);
    
    cs *Operator = cs_spalloc(3, 3, nnz, 1, 0);

    int j;
    for(j = 0; j < nnz; j++){
      Operator->i[j] = I[j];
      Operator->p[j] = J[j];
      Operator->x[j] = val[j];
      Operator-> nz++;
    }

    // Convert a triplet form to compressed column form
    //Operator = cs_compress(Operator);

    // Right hand side
    double b[] = {1, 1, 1, 1};
    
    // Solving Ax = b
    int status = cs_cholsol(Operator, &b[0], 0); // status = 0 means error
    printf("status: %d\n", status);
    
    // use csparse
    
    // print result
  
    return 0;
}

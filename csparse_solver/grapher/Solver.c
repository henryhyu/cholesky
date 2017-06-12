# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include "csparse.h"
# include "mmio.h"
# include "json.h"

void put_json_entry(json_object *jobj, const char *name, const char *y_label, const char *x_label) {
    json_object *entry = json_object_new_object();
    json_object *xlab = json_object_new_string(x_label);
    json_object *ylab = json_object_new_string(y_label);
    json_object *xdata = json_object_new_array();
    json_object *ydata = json_object_new_array();
    json_object *files = json_object_new_array();

    json_object_object_add(entry,"x_label", xlab);
    json_object_object_add(entry,"y_label", ylab);
    json_object_object_add(entry,"x_data", xdata);
    json_object_object_add(entry,"y_data", ydata);
    json_object_object_add(entry,"files", files);
    json_object_object_add(jobj, name, entry);
}

void update_ydata(json_object *jobj, double value, const char *name, const char *fname) {
    json_object *entry = json_object_new_object();
    json_object *ydata = json_object_new_object();
    json_object *files = json_object_new_object();
    
    json_object_object_get_ex(jobj, name, &entry);
    json_object_object_get_ex(entry, "y_data", &ydata);
    json_object_object_get_ex(entry, "files", &files);

    json_object *jvalue = json_object_new_double(value);
    json_object_array_add(ydata, jvalue);

    json_object *file_name = json_object_new_string(fname);
    json_object_array_add(files, file_name);
}

void update_xdata(json_object *jobj, double value, const char *name) {
    json_object *entry = json_object_new_object();
    json_object *xdata = json_object_new_object();
    json_object *files = json_object_new_object();
    
    json_object_object_get_ex(jobj, name, &entry);
    json_object_object_get_ex(entry, "x_data", &xdata);

    json_object *jvalue = json_object_new_double(value);
    json_object_array_add(xdata, jvalue);
}



int main (int argc, char *argv[]) {

    int ret_code;
    FILE *f;
    int row, col, nz, i;  
    int *I, *J, n;
    double *val, *b;
    int order=0;

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
    if (ok) {
        start = clock();
        cs_ipvec (n, S->Pinv, b, x) ;   /* x = P*b */
        cs_lsolve (N->L, x) ;       /* x = L\x */
        cs_ltsolve (N->L, x) ;      /* x = L'\x */
        cs_pvec (n, S->Pinv, x, b) ;    /* b = P'*x */
        solve_time = clock()-start;
    }
    total_time = clock()-start2;

    // below is hard coded*
    json_object *jobj = json_object_from_file(argv[2]);
    if (jobj == NULL) {
        jobj = json_object_new_object();
        put_json_entry(jobj, "Ordering", "Time (s)", "Num Edges");
        put_json_entry(jobj, "Factoring", "Time (s)", "Num Edges");
        put_json_entry(jobj, "TriangularSolve", "Time (s)", "Num Edges");
        put_json_entry(jobj, "TotalTime", "Time (s)", "Num Edges");
        put_json_entry(jobj, "nnz", "Non-zeroes", "Num Edges");
        put_json_entry(jobj, "nnzRatio", "Non-zeroes ratio", "Num Edges");
    }
    update_ydata(jobj, symb_time*1.0/CLOCKS_PER_SEC, "Ordering", argv[1]);
    update_ydata(jobj, num_time*1.0/CLOCKS_PER_SEC, "Factoring", argv[1]);
    update_ydata(jobj, solve_time*1.0/CLOCKS_PER_SEC, "TriangularSolve", argv[1]);
    update_ydata(jobj, total_time*1.0/CLOCKS_PER_SEC, "TotalTime", argv[1]);
    update_ydata(jobj, 1.0*N->L->nzmax/A->nzmax, "nnzRatio", argv[1]);

    update_ydata(jobj, 1.0*A->nzmax, "nnz", argv[1]);
    update_xdata(jobj, 1.0*N->L->nzmax, "nnz");

    json_object_to_file(argv[2], jobj);
    // printf("jobj from str:\n---\n%s\n---\n", json_object_to_json_string_ext(jobj, JSON_C_TO_STRING_SPACED | JSON_C_TO_STRING_PRETTY));
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

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <stdbool.h>
#include <omp.h>
#include "mmio.c"

struct timespec t_start, t_end;


u_int32_t coo2csc(
  u_int32_t       * const row,       /*!< CSC row start indices */
  u_int32_t       * const col,       /*!< CSC column indices */
  u_int32_t const * const row_coo,   /*!< COO row indices */
  u_int32_t const * const col_coo,   /*!< COO column indices */
  u_int32_t const         nnz,       /*!< Number of nonzero elements */
  u_int32_t const         n,         /*!< Number of rows/columns */
  u_int32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

  // ----- cannot assume that input is already 0!
    for(u_int32_t l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
    for(u_int32_t l = 0; l < nnz; l++)
        col[col_coo[l] - isOneBased]++;

  // ----- cumulative sum
    for(u_int32_t i = 0, cumsum = 0; i < n; i++) {
        u_int32_t temp = col[i];
        col[i] = cumsum;
        cumsum += temp;
    }
    col[n] = nnz;
  // ----- copy the row indices to the correct place
    for(u_int32_t l = 0; l < nnz; l++) {
        u_int32_t col_l;
        col_l = col_coo[l] - isOneBased;

        u_int32_t dst = col[col_l];
        row[dst] = row_coo[l] - isOneBased;

        col[col_l]++;
    }
  // ----- revert the column pointers
    for(u_int32_t i = 0, last = 0; i < n; i++) {
        u_int32_t temp = col[i];
        col[i] = last;
        last = temp;
    }

    return n;

}


/* Reads a MMfile */
u_int32_t cooReader(char* name, u_int32_t* I, u_int32_t* J, u_int32_t* II, u_int32_t* JJ){

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   
    int i;
    double *val;


    if ((f = fopen( name, "r")) == NULL) 
        exit(1);
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */

    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d \n", &I[i], &J[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
   // for (i=0; i<nz; i++)
      //  fprintf(stdout, "%d %d \n", I[i]+1, J[i]+1);

	
	//return converter(I,J,val,nz,nz,nz,II,JJ);
	

	printf("nzz=%d\n",nz);
  	return coo2csc(II, JJ, I, J,nz, N,0);
}



int* V4_omp(int *row, int *col, int N, int threadnum){

    omp_set_num_threads(threadnum);
       
    int tr = 0;
    int *c3 = (int *)calloc(N, sizeof(int));

     
    clock_gettime(CLOCK_REALTIME, &t_start);

    #pragma omp parallel
    {
    #pragma omp for nowait schedule(dynamic) reduction(+: tr) 
    for(int i=0; i<N; i++){

        for(int j=row[i]; j<row[i+1]; j++){

            int common = 0;
            int p1     = 0;
            int p2     = 0;

            while(row[i]+p1 < row[i+1] && row[col[j]]+p2 < row[col[j]+1]){

                if(col[row[i]+p1] < col[row[col[j]]+p2]){
                    p1++;
                }
                else if(col[row[i]+p1] > col[row[col[j]]+p2]){
                    p2++;
                }
                else {

                    c3[col[row[i]+p1]]++;
                    c3[i]++;
                    c3[col[j]]++;

                    p1++;
                    p2++;
                    common++;
                }
            }
            tr += common;
        }
    }
    }
    
    clock_gettime(CLOCK_REALTIME, &t_end);

    double duration = ((t_end.tv_sec - t_start.tv_sec) * 1000000 + (t_end.tv_nsec - t_start.tv_nsec) / 1000) / 1000000.0;
    printf("~ Duration : %f sec\n", duration);
    printf("~ Triangles: %d\n", tr);
   
    return c3;
}



int main(int argc, char* argv[]){
    
    char *str           = argv[1];
    int combinationsNum = atoi(argv[2]);
    int rowsNum         = atoi(argv[3]);
    int threadnum       = atoi(argv[4]);

    u_int32_t *I;
    u_int32_t *J;

    int *CSCrows = (int *) malloc(combinationsNum * sizeof(int));
    int *CSCcols = (int *) malloc((rowsNum + 1)   * sizeof(int));

    u_int32_t rowptrSize = cooReader(str, I, J, CSCrows, CSCcols);

    printf("Graph: %s\n", str);
    V4_omp(CSCcols, CSCrows, rowptrSize, threadnum);

    return 0;
}
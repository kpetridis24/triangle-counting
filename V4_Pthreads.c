#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <stdbool.h>
#include <pthread.h>
#include "mmio.c"

struct timespec t_start, t_end;

double TOTAL_TIME = 0;
int TRIANGLES     = 0;


struct thread_args{
    
    int *rows;
    int *cols;
    int rowPtrSize;
    int tid;
    int threadnum;

};



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



void *V4_pthreads(void *arg){
       
    struct thread_args *targs = arg;
    int *row   = targs -> rows;
    int *col   = targs -> cols;
    int N      = targs -> rowPtrSize;
    int id     = targs -> tid;
    int nthrds = targs -> threadnum;

    int *c3 = (int *)calloc(N, sizeof(int));

    int tr          = 0;
    int bound       = 0;
    int extrabound  = 0;
    int threadStart = 0;
    int threadEnd   = 0;
    int t_ptr       = 0;

    if(nthrds > N)
    {
        printf("Number of threads greater than number of iterations!");
        exit(0);
    }


    /************************************
    * DISTRIBUTE THE ITERATIONS EQUALLY *
    ************************************/

   if(N % nthrds == 0)
    {
        bound      = N / nthrds;
        extrabound = 0;
    }
    else
    {
        bound      = N / nthrds;
        extrabound = N % nthrds;
    }

    if(t_ptr != nthrds - 1 || nthrds == 1)
    {
        threadStart = id * bound;
        threadEnd   = threadStart + bound;
    }

    if(t_ptr == nthrds - 1 && nthrds != 1)
    {
        threadStart = N - extrabound;
        threadEnd   = N;
    }

    t_ptr++;


    
    clock_gettime(CLOCK_REALTIME, &t_start);

    for(int i = threadStart; i < threadEnd; i++){

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
    
    clock_gettime(CLOCK_REALTIME, &t_end);

    double duration = ((t_end.tv_sec - t_start.tv_sec) * 1000000 + (t_end.tv_nsec - t_start.tv_nsec) / 1000) / 1000000.0;
    printf("Thread(%d) duration: %f\tseconds\n", id, duration);
    printf("Thread(%d) detected: %d\ttriangles\n\n", id, tr);

    TOTAL_TIME += duration;
    TRIANGLES  += tr;

    pthread_exit(NULL);

    return (void *)c3;

}



void main(int argc, char *argv[]){
    
    char* str = argv[1];
    int combinationsNum = atoi(argv[2]);
    int rowsNum = atoi(argv[3]);
    int tnum = atoi(argv[4]);

    pthread_t *threads        = (pthread_t *)malloc(tnum * sizeof(pthread_t));
    struct thread_args *targs = (struct thread_args *)malloc(tnum * sizeof(struct thread_args));


    u_int32_t *I;
    u_int32_t *J;

    u_int32_t *CSCrows = (u_int32_t *) malloc(combinationsNum * sizeof(u_int32_t));
    u_int32_t *CSCcols = (u_int32_t *) malloc((rowsNum + 1)* sizeof(u_int32_t));

    u_int32_t N  = cooReader(str, I, J, CSCrows, CSCcols);
    

    /* Initialization */

    for(int i=0, k=0; i<tnum; i++)
    {
        targs[i].rows       = CSCcols;
        targs[i].cols       = CSCrows;
        targs[i].tid        = i;
        targs[i].rowPtrSize = N;
        targs[i].threadnum  = tnum;

        pthread_create(&threads[i], NULL, V4_pthreads, (void *)&targs[i]);
    }
    
    for(int l=0; l<tnum; l++)
        pthread_join(threads[l], NULL);


    printf("Total duration : %lf\n", TOTAL_TIME);
    printf("Total triangles: %d\n", TRIANGLES);
   
}
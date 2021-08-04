#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "gemm.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*gemm.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int ni, int nj, int nk, double *alpha, double *beta, double C[1000][1100], double A[1000][1200], double B[1200][1100]) {
   int i, j;
   *alpha = 1.5;
   *beta = 1.2;
   for(i = 0; i < ni; i++)
      for(j = 0; j < nj; j++)
         C[i][j] = (double) ((i * j + 1) % ni) / ni;
   for(i = 0; i < ni; i++)
      for(j = 0; j < nk; j++)
         A[i][j] = (double) (i * (j + 1) % nk) / nk;
   for(i = 0; i < nk; i++)
      for(j = 0; j < nj; j++)
         B[i][j] = (double) (i * (j + 2) % nj) / nj;
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int ni, int nj, double C[1000][1100]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "C");
   for(i = 0; i < ni; i++)
      for(j = 0; j < nj; j++) {
         if((i * ni + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", C[i][j]);
      }
   fprintf(stderr, "\nend   dump: %s\n", "C");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_gemm(int ni, int nj, int nk, double alpha, double beta, double C[1000][1100], double A[1000][1200], double B[1200][1100]) {
   int i, j, k;
   #pragma omp parallel for default(shared) private(i, j, k) firstprivate(ni, nj, beta, nk, alpha, A, B)
   for(i = 0; i < ni; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(nj, i, beta)
      for(j = 0; j < nj; j++)
         C[i][j] *= beta;
      // #pragma omp parallel for default(shared) private(k, j) firstprivate(nk, nj, alpha, i, A, B)
      for(k = 0; k < nk; k++) {
         // #pragma omp parallel for default(shared) private(j) firstprivate(nj, alpha, i, k, A, B)
         for(j = 0; j < nj; j++)
            C[i][j] += alpha * A[i][k] * B[k][j];
      }
   }
}

int main(int argc, char **argv) {
   /*Retrieve problem size.*/
   int ni = 1000;
   int nj = 1100;
   int nk = 1200;
   /*Variable declaration/allocation.*/
   double alpha;
   double beta;
   double (*C)[1000][1100];
   C = (double (*)[1000][1100]) polybench_alloc_data((1000 + 0) * (1100 + 0), sizeof(double));
   ;
   double (*A)[1000][1200];
   A = (double (*)[1000][1200]) polybench_alloc_data((1000 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*B)[1200][1100];
   B = (double (*)[1200][1100]) polybench_alloc_data((1200 + 0) * (1100 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(ni, nj, nk, &alpha, &beta, *C, *A, *B);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_gemm(ni, nj, nk, alpha, beta, *C, *A, *B);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(ni, nj, *C);
   /*Be clean.*/
   free((void *) C);
   ;
   free((void *) A);
   ;
   free((void *) B);
   ;
   
   return 0;
}

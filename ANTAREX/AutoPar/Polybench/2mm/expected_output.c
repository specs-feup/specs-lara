#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "2mm.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*2mm.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int ni, int nj, int nk, int nl, double *alpha, double *beta, double A[800][1100], double B[1100][900], double C[900][1200], double D[800][1200]) {
   int i, j;
   *alpha = 1.5;
   *beta = 1.2;
   for(i = 0; i < ni; i++)
      for(j = 0; j < nk; j++)
         A[i][j] = (double) ((i * j + 1) % ni) / ni;
   for(i = 0; i < nk; i++)
      for(j = 0; j < nj; j++)
         B[i][j] = (double) (i * (j + 1) % nj) / nj;
   for(i = 0; i < nj; i++)
      for(j = 0; j < nl; j++)
         C[i][j] = (double) ((i * (j + 3) + 1) % nl) / nl;
   for(i = 0; i < ni; i++)
      for(j = 0; j < nl; j++)
         D[i][j] = (double) (i * (j + 2) % nk) / nk;
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int ni, int nl, double D[800][1200]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "D");
   for(i = 0; i < ni; i++)
      for(j = 0; j < nl; j++) {
         if((i * ni + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", D[i][j]);
      }
   fprintf(stderr, "\nend   dump: %s\n", "D");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_2mm(int ni, int nj, int nk, int nl, double alpha, double beta, double tmp[800][900], double A[800][1100], double B[1100][900], double C[900][1200], double D[800][1200]) {
   int i, j, k;
   #pragma omp parallel for default(shared) private(i, j, k) firstprivate(ni, nj, nk, alpha, A, B)
   for(i = 0; i < ni; i++) {
      // #pragma omp parallel for default(shared) private(j, k) firstprivate(nj, i, nk, alpha, A, B)
      for(j = 0; j < nj; j++) {
         tmp[i][j] = 0.0;
         // #pragma omp parallel for default(shared) private(k) firstprivate(nk, alpha, i, j, A, B) reduction(+ : tmp[i][j])
         for(k = 0; k < nk; ++k)
            tmp[i][j] += alpha * A[i][k] * B[k][j];
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k) firstprivate(ni, nl, beta, nj, tmp, C)
   for(i = 0; i < ni; i++) {
      // #pragma omp parallel for default(shared) private(j, k) firstprivate(nl, i, beta, nj, tmp, C)
      for(j = 0; j < nl; j++) {
         D[i][j] *= beta;
         // #pragma omp parallel for default(shared) private(k) firstprivate(nj, i, j, tmp, C) reduction(+ : D[i][j])
         for(k = 0; k < nj; ++k)
            D[i][j] += tmp[i][k] * C[k][j];
      }
   }
}

int main(int argc, char **argv) {
   /*Retrieve problem size.*/
   int ni = 800;
   int nj = 900;
   int nk = 1100;
   int nl = 1200;
   /*Variable declaration/allocation.*/
   double alpha;
   double beta;
   double (*tmp)[800][900];
   tmp = (double (*)[800][900]) polybench_alloc_data((800 + 0) * (900 + 0), sizeof(double));
   ;
   double (*A)[800][1100];
   A = (double (*)[800][1100]) polybench_alloc_data((800 + 0) * (1100 + 0), sizeof(double));
   ;
   double (*B)[1100][900];
   B = (double (*)[1100][900]) polybench_alloc_data((1100 + 0) * (900 + 0), sizeof(double));
   ;
   double (*C)[900][1200];
   C = (double (*)[900][1200]) polybench_alloc_data((900 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*D)[800][1200];
   D = (double (*)[800][1200]) polybench_alloc_data((800 + 0) * (1200 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(ni, nj, nk, nl, &alpha, &beta, *A, *B, *C, *D);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_2mm(ni, nj, nk, nl, alpha, beta, *tmp, *A, *B, *C, *D);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(ni, nl, *D);
   /*Be clean.*/
   free((void *) tmp);
   ;
   free((void *) A);
   ;
   free((void *) B);
   ;
   free((void *) C);
   ;
   free((void *) D);
   ;
   
   return 0;
}

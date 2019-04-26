#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "syr2k.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*syr2k.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int n, int m, double *alpha, double *beta, double C[1200][1200], double A[1200][1000], double B[1200][1000]) {
   int i, j;
   *alpha = 1.5;
   *beta = 1.2;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n, m)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(m, i, n)
      for(j = 0; j < m; j++) {
         A[i][j] = (double) ((i * j + 1) % n) / n;
         B[i][j] = (double) ((i * j + 2) % m) / m;
      }
   }
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n, m)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i, m)
      for(j = 0; j < n; j++) {
         C[i][j] = (double) ((i * j + 3) % n) / m;
      }
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int n, double C[1200][1200]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "C");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#59{fprintf(stderr, "\n")}
   fprintf#61{fprintf(stderr, "%0.2lf ", C[i][j])}
   ****************************************/
   for(i = 0; i < n; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#59{fprintf(stderr, "\n")}
      fprintf#61{fprintf(stderr, "%0.2lf ", C[i][j])}
      ****************************************/
      for(j = 0; j < n; j++) {
         if((i * n + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", C[i][j]);
      }
   }
   fprintf(stderr, "\nend   dump: %s\n", "C");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_syr2k(int n, int m, double alpha, double beta, double C[1200][1200], double A[1200][1000], double B[1200][1000]) {
   int i, j, k;
   #pragma omp parallel for default(shared) private(i, j, k) firstprivate(n, beta, m, alpha)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(i, beta)
      for(j = 0; j <= i; j++) C[i][j] *= beta;
      // #pragma omp parallel for default(shared) private(k, j) firstprivate(m, i, alpha)
      for(k = 0; k < m; k++) {
         // #pragma omp parallel for default(shared) private(j) firstprivate(i, k, alpha)
         for(j = 0; j <= i; j++) {
            C[i][j] += A[j][k] * alpha * B[i][k] + B[j][k] * alpha * A[i][k];
         }
      }
   }
}

int main(int argc, char **argv) {
   /*Retrieve problem size.*/
   int n = 1200;
   int m = 1000;
   /*Variable declaration/allocation.*/
   double alpha;
   double beta;
   double (*C)[1200][1200];
   C = (double (*)[1200][1200]) polybench_alloc_data((1200 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*A)[1200][1000];
   A = (double (*)[1200][1000]) polybench_alloc_data((1200 + 0) * (1000 + 0), sizeof(double));
   ;
   double (*B)[1200][1000];
   B = (double (*)[1200][1000]) polybench_alloc_data((1200 + 0) * (1000 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(n, m, &alpha, &beta, *C, *A, *B);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_syr2k(n, m, alpha, beta, *C, *A, *B);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(n, *C);
   /*Be clean.*/
   free((void *) C);
   ;
   free((void *) A);
   ;
   free((void *) B);
   ;
   
   return 0;
}

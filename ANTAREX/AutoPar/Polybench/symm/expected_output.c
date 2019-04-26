#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "symm.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*symm.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int m, int n, double *alpha, double *beta, double C[1000][1200], double A[1000][1000], double B[1000][1200]) {
   int i, j;
   *alpha = 1.5;
   *beta = 1.2;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(m, n)
   for(i = 0; i < m; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i, m)
      for(j = 0; j < n; j++) {
         C[i][j] = (double) ((i + j) % 100) / m;
         B[i][j] = (double) ((n + i - j) % 100) / m;
      }
   }
   #pragma omp parallel for default(shared) private(i, j) firstprivate(m)
   for(i = 0; i < m; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(i, m)
      for(j = 0; j <= i; j++) A[i][j] = (double) ((i + j) % 100) / m;
      // #pragma omp parallel for default(shared) private(j) firstprivate(i, m)
      for(j = i + 1; j < m; j++) A[i][j] = -999; //regions of arrays that should not be used
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int m, int n, double C[1000][1200]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "C");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#64{fprintf(stderr, "\n")}
   fprintf#66{fprintf(stderr, "%0.2lf ", C[i][j])}
   ****************************************/
   for(i = 0; i < m; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#64{fprintf(stderr, "\n")}
      fprintf#66{fprintf(stderr, "%0.2lf ", C[i][j])}
      ****************************************/
      for(j = 0; j < n; j++) {
         if((i * m + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", C[i][j]);
      }
   }
   fprintf(stderr, "\nend   dump: %s\n", "C");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_symm(int m, int n, double alpha, double beta, double C[1000][1200], double A[1000][1000], double B[1000][1200]) {
   int i, j, k;
   double temp2;
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess C	 use : RW
   ****************************************/
   for(i = 0; i < m; i++) {
      #pragma omp parallel for default(shared) private(j, k, temp2) firstprivate(n, i, alpha, beta)
      for(j = 0; j < n; j++) {
         temp2 = 0;
         // #pragma omp parallel for default(shared) private(k) firstprivate(i, alpha, j) reduction(+ : temp2)
         for(k = 0; k < i; k++) {
            C[k][j] += alpha * B[i][j] * A[i][k];
            temp2 += B[k][j] * A[i][k];
         }
         C[i][j] = beta * C[i][j] + alpha * B[i][j] * A[i][i] + alpha * temp2;
      }
   }
}

int main(int argc, char **argv) {
   /*Retrieve problem size.*/
   int m = 1000;
   int n = 1200;
   /*Variable declaration/allocation.*/
   double alpha;
   double beta;
   double (*C)[1000][1200];
   C = (double (*)[1000][1200]) polybench_alloc_data((1000 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*A)[1000][1000];
   A = (double (*)[1000][1000]) polybench_alloc_data((1000 + 0) * (1000 + 0), sizeof(double));
   ;
   double (*B)[1000][1200];
   B = (double (*)[1000][1200]) polybench_alloc_data((1000 + 0) * (1200 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(m, n, &alpha, &beta, *C, *A, *B);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_symm(m, n, alpha, beta, *C, *A, *B);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(m, n, *C);
   /*Be clean.*/
   free((void *) C);
   ;
   free((void *) A);
   ;
   free((void *) B);
   ;
   
   return 0;
}

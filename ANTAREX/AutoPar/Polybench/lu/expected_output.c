#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "lu.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*lu.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int n, double A[2000][2000]) {
   int i, j;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(i, n)
      for(j = 0; j <= i; j++) A[i][j] = (double) (-j % n) / n + 1;
      // #pragma omp parallel for default(shared) private(j) firstprivate(i, n)
      for(j = i + 1; j < n; j++) {
         A[i][j] = 0;
      }
      A[i][i] = 1;
   }
   /*Make the matrix positive semi-definite.*/
   /*not necessary for LU, but using same code as cholesky*/
   int r, s, t;
   double (*B)[2000][2000];
   B = (double (*)[2000][2000]) polybench_alloc_data((2000 + 0) * (2000 + 0), sizeof(double));
   ;
   #pragma omp parallel for default(shared) private(r, s) firstprivate(n)
   for(r = 0; r < n; ++r) {
      // #pragma omp parallel for default(shared) private(s) firstprivate(n, r)
      for(s = 0; s < n; ++s) (*B)[r][s] = 0;
   }
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess (*B)	 use : RW
   ****************************************/
   for(t = 0; t < n; ++t) {
      #pragma omp parallel for default(shared) private(r, s) firstprivate(n, t)
      for(r = 0; r < n; ++r) {
         // #pragma omp parallel for default(shared) private(s) firstprivate(n, r, t)
         for(s = 0; s < n; ++s) (*B)[r][s] += A[r][t] * A[s][t];
      }
   }
   #pragma omp parallel for default(shared) private(r, s) firstprivate(n)
   for(r = 0; r < n; ++r) {
      // #pragma omp parallel for default(shared) private(s) firstprivate(n, r)
      for(s = 0; s < n; ++s) A[r][s] = (*B)[r][s];
   }
   free((void *) B);
   ;
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int n, double A[2000][2000]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "A");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#92{fprintf(stderr, "\n")}
   fprintf#94{fprintf(stderr, "%0.2lf ", A[i][j])}
   ****************************************/
   for(i = 0; i < n; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#92{fprintf(stderr, "\n")}
      fprintf#94{fprintf(stderr, "%0.2lf ", A[i][j])}
      ****************************************/
      for(j = 0; j < n; j++) {
         if((i * n + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", A[i][j]);
      }
   }
   fprintf(stderr, "\nend   dump: %s\n", "A");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_lu(int n, double A[2000][2000]) {
   int i, j, k;
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess A	 use : RW
   ****************************************/
   for(i = 0; i < n; i++) {
      /*************** Clava msgError **************
      unsolved dependency for arrayAccess A	 use : RWR
      ****************************************/
      for(j = 0; j < i; j++) {
         /*************** Clava msgError **************
         unsolved dependency for arrayAccess A	 use : RW
         ****************************************/
         for(k = 0; k < j; k++) {
            A[i][j] -= A[i][k] * A[k][j];
         }
         A[i][j] /= A[j][j];
      }
      /*************** Clava msgError **************
      unsolved dependency for arrayAccess A	 use : RW
      ****************************************/
      for(j = i; j < n; j++) {
         /*************** Clava msgError **************
         unsolved dependency for arrayAccess A	 use : RW
         ****************************************/
         for(k = 0; k < i; k++) {
            A[i][j] -= A[i][k] * A[k][j];
         }
      }
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int n = 2000;
   /*Variable declaration/allocation.*/
   double (*A)[2000][2000];
   A = (double (*)[2000][2000]) polybench_alloc_data((2000 + 0) * (2000 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(n, *A);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_lu(n, *A);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(n, *A);
   /*Be clean.*/
   free((void *) A);
   ;
   
   return 0;
}

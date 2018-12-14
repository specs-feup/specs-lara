#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "ludcmp.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*ludcmp.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int n, double A[2000][2000], double b[2000], double x[2000], double y[2000]) {
   int i, j;
   double fn = (double) n;
   #pragma omp parallel for default(shared) private(i) firstprivate(n, fn)
   for(i = 0; i < n; i++) {
      x[i] = 0;
      y[i] = 0;
      b[i] = (i + 1) / fn / 2.0 + 4;
   }
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
static void print_array(int n, double x[2000]) {
   int i;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "x");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#97{fprintf(stderr, "\n")}
   fprintf#99{fprintf(stderr, "%0.2lf ", x[i])}
   ****************************************/
   for(i = 0; i < n; i++) {
      if(i % 20 == 0) fprintf(stderr, "\n");
      fprintf(stderr, "%0.2lf ", x[i]);
   }
   fprintf(stderr, "\nend   dump: %s\n", "x");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_ludcmp(int n, double A[2000][2000], double b[2000], double x[2000], double y[2000]) {
   int i, j, k;
   double w;
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess A	 use : RW
   ****************************************/
   for(i = 0; i < n; i++) {
      /*************** Clava msgError **************
      unsolved dependency for arrayAccess A	 use : RW
      ****************************************/
      for(j = 0; j < i; j++) {
         w = A[i][j];
         #pragma omp parallel for default(shared) private(k) firstprivate(j, i) reduction(- : w)
         for(k = 0; k < j; k++) {
            w -= A[i][k] * A[k][j];
         }
         A[i][j] = w / A[j][j];
      }
      /*************** Clava msgError **************
      unsolved dependency for arrayAccess A	 use : RW
      ****************************************/
      for(j = i; j < n; j++) {
         w = A[i][j];
         #pragma omp parallel for default(shared) private(k) firstprivate(i, j) reduction(- : w)
         for(k = 0; k < i; k++) {
            w -= A[i][k] * A[k][j];
         }
         A[i][j] = w;
      }
   }
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess y	 use : RW
   ****************************************/
   for(i = 0; i < n; i++) {
      w = b[i];
      #pragma omp parallel for default(shared) private(j) firstprivate(i) reduction(- : w)
      for(j = 0; j < i; j++) w -= A[i][j] * y[j];
      y[i] = w;
   }
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess x	 use : RW
   ****************************************/
   for(i = n - 1; i >= 0; i--) {
      w = y[i];
      #pragma omp parallel for default(shared) private(j) firstprivate(i, n) reduction(- : w)
      for(j = i + 1; j < n; j++) w -= A[i][j] * x[j];
      x[i] = w / A[i][i];
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int n = 2000;
   double (*A)[2000][2000];
   /*Variable declaration/allocation.*/
   A = (double (*)[2000][2000]) polybench_alloc_data((2000 + 0) * (2000 + 0), sizeof(double));
   ;
   double (*b)[2000];
   b = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*x)[2000];
   x = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*y)[2000];
   y = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(n, *A, *b, *x, *y);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_ludcmp(n, *A, *b, *x, *y);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(n, *x);
   free((void *) A);
   /*Be clean.*/
   ;
   free((void *) b);
   ;
   free((void *) x);
   ;
   free((void *) y);
   ;
   
   return 0;
}

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "gemver.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*gemver.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int n, double * alpha, double * beta, double A[2000][2000], double u1[2000], double v1[2000], double u2[2000], double v2[2000], double w[2000], double x[2000], double y[2000], double z[2000]) {
   int i, j;
   *alpha = 1.5;
   *beta = 1.2;
   double fn = (double) n;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n, fn)
   for(i = 0; i < n; i++) {
      u1[i] = i;
      u2[i] = ((i + 1) / fn) / 2.0;
      v1[i] = ((i + 1) / fn) / 4.0;
      v2[i] = ((i + 1) / fn) / 6.0;
      y[i] = ((i + 1) / fn) / 8.0;
      z[i] = ((i + 1) / fn) / 9.0;
      x[i] = 0.0;
      w[i] = 0.0;
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i)
      for(j = 0; j < n; j++) A[i][j] = (double) (i * j % n) / n;
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int n, double w[2000]) {
   int i;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "w");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#55{fprintf(stderr, "\n")}
   fprintf#57{fprintf(stderr, "%0.2lf ", w[i])}
   ****************************************/
   for(i = 0; i < n; i++) {
      if(i % 20 == 0) fprintf(stderr, "\n");
      fprintf(stderr, "%0.2lf ", w[i]);
   }
   fprintf(stderr, "\nend   dump: %s\n", "w");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_gemver(int n, double alpha, double beta, double A[2000][2000], double u1[2000], double v1[2000], double u2[2000], double v2[2000], double w[2000], double x[2000], double y[2000], double z[2000]) {
   int i, j;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i)
      for(j = 0; j < n; j++) A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];
   }
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n, beta)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, beta, i) reduction(+ : x[i])
      for(j = 0; j < n; j++) x[i] = x[i] + beta * A[j][i] * y[j];
   }
   #pragma omp parallel for default(shared) private(i) firstprivate(n)
   for(i = 0; i < n; i++) x[i] = x[i] + z[i];
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n, alpha)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, alpha, i) reduction(+ : w[i])
      for(j = 0; j < n; j++) w[i] = w[i] + alpha * A[i][j] * x[j];
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int n = 2000;
   /*Variable declaration/allocation.*/
   double alpha;
   double beta;
   double (*A)[2000][2000];
   A = (double (*)[2000][2000]) polybench_alloc_data((2000 + 0) * (2000 + 0), sizeof(double));
   ;
   double (*u1)[2000];
   u1 = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*v1)[2000];
   v1 = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*u2)[2000];
   u2 = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*v2)[2000];
   v2 = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*w)[2000];
   w = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*x)[2000];
   x = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*y)[2000];
   y = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*z)[2000];
   z = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(n, &alpha, &beta, *A, *u1, *v1, *u2, *v2, *w, *x, *y, *z);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_gemver(n, alpha, beta, *A, *u1, *v1, *u2, *v2, *w, *x, *y, *z);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(n, *w);
   free((void *) A);
   /*Be clean.*/
   ;
   free((void *) u1);
   ;
   free((void *) v1);
   ;
   free((void *) u2);
   ;
   free((void *) v2);
   ;
   free((void *) w);
   ;
   free((void *) x);
   ;
   free((void *) y);
   ;
   free((void *) z);
   ;
   
   return 0;
}

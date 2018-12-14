#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "atax.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*atax.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int m, int n, double A[1900][2100], double x[2100]) {
   int i, j;
   double fn;
   fn = (double) n;
   #pragma omp parallel for default(shared) private(i) firstprivate(n, fn)
   for(i = 0; i < n; i++) x[i] = 1 + (i / fn);
   #pragma omp parallel for default(shared) private(i, j) firstprivate(m, n)
   for(i = 0; i < m; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i, m)
      for(j = 0; j < n; j++) A[i][j] = (double) ((i + j) % n) / (5 * m);
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int n, double y[2100]) {
   int i;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "y");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#51{fprintf(stderr, "\n")}
   fprintf#53{fprintf(stderr, "%0.2lf ", y[i])}
   ****************************************/
   for(i = 0; i < n; i++) {
      if(i % 20 == 0) fprintf(stderr, "\n");
      fprintf(stderr, "%0.2lf ", y[i]);
   }
   fprintf(stderr, "\nend   dump: %s\n", "y");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_atax(int m, int n, double A[1900][2100], double x[2100], double y[2100], double tmp[1900]) {
   int i, j;
   #pragma omp parallel for default(shared) private(i) firstprivate(n)
   for(i = 0; i < n; i++) y[i] = 0;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(m, n) reduction(+ : y[:2100])
   for(i = 0; i < m; i++) {
      tmp[i] = 0.0;
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i) reduction(+ : tmp[i])
      for(j = 0; j < n; j++) tmp[i] = tmp[i] + A[i][j] * x[j];
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i)
      for(j = 0; j < n; j++) y[j] = y[j] + A[i][j] * tmp[i];
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int m = 1900;
   int n = 2100;
   double (*A)[1900][2100];
   /*Variable declaration/allocation.*/
   A = (double (*)[1900][2100]) polybench_alloc_data((1900 + 0) * (2100 + 0), sizeof(double));
   ;
   double (*x)[2100];
   x = (double (*)[2100]) polybench_alloc_data(2100 + 0, sizeof(double));
   ;
   double (*y)[2100];
   y = (double (*)[2100]) polybench_alloc_data(2100 + 0, sizeof(double));
   ;
   double (*tmp)[1900];
   tmp = (double (*)[1900]) polybench_alloc_data(1900 + 0, sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(m, n, *A, *x);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_atax(m, n, *A, *x, *y, *tmp);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(n, *y);
   free((void *) A);
   /*Be clean.*/
   ;
   free((void *) x);
   ;
   free((void *) y);
   ;
   free((void *) tmp);
   ;
   
   return 0;
}

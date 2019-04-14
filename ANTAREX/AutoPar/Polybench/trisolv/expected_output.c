#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "trisolv.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*trisolv.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int n, double L[2000][2000], double x[2000], double b[2000]) {
   int i, j;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n)
   for(i = 0; i < n; i++) {
      x[i] = -999;
      b[i] = i;
      // #pragma omp parallel for default(shared) private(j) firstprivate(i, n)
      for(j = 0; j <= i; j++) L[i][j] = (double) (i + n - j + 1) * 2 / n;
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int n, double x[2000]) {
   int i;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "x");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#45{fprintf(stderr, "%0.2lf ", x[i])}
   fprintf#47{fprintf(stderr, "\n")}
   ****************************************/
   for(i = 0; i < n; i++) {
      fprintf(stderr, "%0.2lf ", x[i]);
      if(i % 20 == 0) fprintf(stderr, "\n");
   }
   fprintf(stderr, "\nend   dump: %s\n", "x");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_trisolv(int n, double L[2000][2000], double x[2000], double b[2000]) {
   int i, j;
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess x	 use : W
   ****************************************/
   for(i = 0; i < n; i++) {
      x[i] = b[i];
      /*************** Clava msgError **************
      unsolved dependency for arrayAccess x	 use : RW
      ****************************************/
      for(j = 0; j < i; j++) x[i] -= L[i][j] * x[j];
      x[i] = x[i] / L[i][i];
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int n = 2000;
   double (*L)[2000][2000];
   /*Variable declaration/allocation.*/
   L = (double (*)[2000][2000]) polybench_alloc_data((2000 + 0) * (2000 + 0), sizeof(double));
   ;
   double (*x)[2000];
   x = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*b)[2000];
   b = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(n, *L, *x, *b);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_trisolv(n, *L, *x, *b);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(n, *x);
   /*Be clean.*/
   free((void *) L);
   ;
   free((void *) x);
   ;
   free((void *) b);
   ;
   
   return 0;
}

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "seidel-2d.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*seidel-2d.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int n, double A[2000][2000]) {
   int i, j;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i)
      for(j = 0; j < n; j++) A[i][j] = ((double) i * (j + 2) + 2) / n;
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int n, double A[2000][2000]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "A");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#47{fprintf(stderr, "\n")}
   fprintf#49{fprintf(stderr, "%0.2lf ", A[i][j])}
   ****************************************/
   for(i = 0; i < n; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#47{fprintf(stderr, "\n")}
      fprintf#49{fprintf(stderr, "%0.2lf ", A[i][j])}
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
static void kernel_seidel_2d(int tsteps, int n, double A[2000][2000]) {
   int t, i, j;
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess A	 use : RW
   ****************************************/
   for(t = 0; t <= tsteps - 1; t++) {
      /*************** Clava msgError **************
      unsolved dependency for arrayAccess A	 use : RW
      ****************************************/
      for(i = 1; i <= n - 2; i++) {
         /*************** Clava msgError **************
         unsolved dependency for arrayAccess A	 use : RW
         ****************************************/
         for(j = 1; j <= n - 2; j++) A[i][j] = (A[i - 1][j - 1] + A[i - 1][j] + A[i - 1][j + 1] + A[i][j - 1] + A[i][j] + A[i][j + 1] + A[i + 1][j - 1] + A[i + 1][j] + A[i + 1][j + 1]) / 9.0;
      }
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int n = 2000;
   int tsteps = 500;
   double (*A)[2000][2000];
   /*Variable declaration/allocation.*/
   A = (double (*)[2000][2000]) polybench_alloc_data((2000 + 0) * (2000 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(n, *A);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_seidel_2d(tsteps, n, *A);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(n, *A);
   free((void *) A);
   /*Be clean.*/
   ;
   
   return 0;
}

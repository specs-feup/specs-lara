#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "heat-3d.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*heat-3d.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int n, double A[120][120][120], double B[120][120][120]) {
   int i, j, k;
   #pragma omp parallel for default(shared) private(i, j, k) firstprivate(n)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j, k) firstprivate(n, i)
      for(j = 0; j < n; j++) {
         // #pragma omp parallel for default(shared) private(k) firstprivate(n, i, j)
         for(k = 0; k < n; k++) A[i][j][k] = B[i][j][k] = (double) (i + j + (n - k)) * 10 / (n);
      }
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int n, double A[120][120][120]) {
   int i, j, k;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "A");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#54{fprintf(stderr, "\n")}
   fprintf#56{fprintf(stderr, "%0.2lf ", A[i][j][k])}
   ****************************************/
   for(i = 0; i < n; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#54{fprintf(stderr, "\n")}
      fprintf#56{fprintf(stderr, "%0.2lf ", A[i][j][k])}
      ****************************************/
      for(j = 0; j < n; j++) {
         /*************** Clava msgError **************
         Variables Access as passed arguments Can not be traced inside of function calls :
         fprintf#54{fprintf(stderr, "\n")}
         fprintf#56{fprintf(stderr, "%0.2lf ", A[i][j][k])}
         ****************************************/
         for(k = 0; k < n; k++) {
            if((i * n * n + j * n + k) % 20 == 0) fprintf(stderr, "\n");
            fprintf(stderr, "%0.2lf ", A[i][j][k]);
         }
      }
   }
   fprintf(stderr, "\nend   dump: %s\n", "A");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_heat_3d(int tsteps, int n, double A[120][120][120], double B[120][120][120]) {
   int t, i, j, k;
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess A	 use : RW
   ****************************************/
   for(t = 1; t <= 500; t++) {
      #pragma omp parallel for default(shared) private(i, j, k) firstprivate(n)
      for(i = 1; i < n - 1; i++) {
         // #pragma omp parallel for default(shared) private(j, k) firstprivate(n, i)
         for(j = 1; j < n - 1; j++) {
            // #pragma omp parallel for default(shared) private(k) firstprivate(n, i, j)
            for(k = 1; k < n - 1; k++) {
               B[i][j][k] = 0.125 * (A[i + 1][j][k] - 2.0 * A[i][j][k] + A[i - 1][j][k]) + 0.125 * (A[i][j + 1][k] - 2.0 * A[i][j][k] + A[i][j - 1][k]) + 0.125 * (A[i][j][k + 1] - 2.0 * A[i][j][k] + A[i][j][k - 1]) + A[i][j][k];
            }
         }
      }
      #pragma omp parallel for default(shared) private(i, j, k) firstprivate(n)
      for(i = 1; i < n - 1; i++) {
         // #pragma omp parallel for default(shared) private(j, k) firstprivate(n, i)
         for(j = 1; j < n - 1; j++) {
            // #pragma omp parallel for default(shared) private(k) firstprivate(n, i, j)
            for(k = 1; k < n - 1; k++) {
               A[i][j][k] = 0.125 * (B[i + 1][j][k] - 2.0 * B[i][j][k] + B[i - 1][j][k]) + 0.125 * (B[i][j + 1][k] - 2.0 * B[i][j][k] + B[i][j - 1][k]) + 0.125 * (B[i][j][k + 1] - 2.0 * B[i][j][k] + B[i][j][k - 1]) + B[i][j][k];
            }
         }
      }
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int n = 120;
   int tsteps = 500;
   double (*A)[120][120][120];
   /*Variable declaration/allocation.*/
   A = (double (*)[120][120][120]) polybench_alloc_data((120 + 0) * (120 + 0) * (120 + 0), sizeof(double));
   ;
   double (*B)[120][120][120];
   B = (double (*)[120][120][120]) polybench_alloc_data((120 + 0) * (120 + 0) * (120 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(n, *A, *B);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_heat_3d(tsteps, n, *A, *B);
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

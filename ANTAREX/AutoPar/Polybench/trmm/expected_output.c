#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "trmm.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*trmm.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int m, int n, double *alpha, double A[1000][1000], double B[1000][1200]) {
   int i, j;
   *alpha = 1.5;
   for(i = 0; i < m; i++) {
      for(j = 0; j < i; j++) {
         A[i][j] = (double) ((i + j) % m) / m;
      }
      A[i][i] = 1.0;
      for(j = 0; j < n; j++) {
         B[i][j] = (double) ((n + (i - j)) % n) / n;
      }
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int m, int n, double B[1000][1200]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "B");
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++) {
         if((i * m + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", B[i][j]);
      }
   fprintf(stderr, "\nend   dump: %s\n", "B");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_trmm(int m, int n, double alpha, double A[1000][1000], double B[1000][1200]) {
   int i, j, k;
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess B	 use : RW
   ****************************************/
   for(i = 0; i < m; i++) {
      #pragma omp parallel for default(shared) private(j, k) firstprivate(n, i, m, alpha, A)
      for(j = 0; j < n; j++) {
         /*************** Clava msgError **************
         unsolved dependency for arrayAccess B	 use : RW
         ****************************************/
         for(k = i + 1; k < m; k++)
            B[i][j] += A[k][i] * B[k][j];
         B[i][j] = alpha * B[i][j];
      }
   }
}

int main(int argc, char **argv) {
   /*Retrieve problem size.*/
   int m = 1000;
   int n = 1200;
   /*Variable declaration/allocation.*/
   double alpha;
   double (*A)[1000][1000];
   A = (double (*)[1000][1000]) polybench_alloc_data((1000 + 0) * (1000 + 0), sizeof(double));
   ;
   double (*B)[1000][1200];
   B = (double (*)[1000][1200]) polybench_alloc_data((1000 + 0) * (1200 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(m, n, &alpha, *A, *B);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_trmm(m, n, alpha, *A, *B);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(m, n, *B);
   /*Be clean.*/
   free((void *) A);
   ;
   free((void *) B);
   ;
   
   return 0;
}

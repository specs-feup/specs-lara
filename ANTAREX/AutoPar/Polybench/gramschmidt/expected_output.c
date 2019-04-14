#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "gramschmidt.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*gramschmidt.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int m, int n, double A[1000][1200], double R[1200][1200], double Q[1000][1200]) {
   int i, j;
   #pragma omp parallel for default(shared) private(i, j) firstprivate(m, n)
   for(i = 0; i < m; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i, m)
      for(j = 0; j < n; j++) {
         A[i][j] = (((double) ((i * j) % m) / m) * 100) + 10;
         Q[i][j] = 0.0;
      }
   }
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(n, i)
      for(j = 0; j < n; j++) R[i][j] = 0.0;
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int m, int n, double A[1000][1200], double R[1200][1200], double Q[1000][1200]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "R");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#57{fprintf(stderr, "\n")}
   fprintf#59{fprintf(stderr, "%0.2lf ", R[i][j])}
   ****************************************/
   for(i = 0; i < n; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#57{fprintf(stderr, "\n")}
      fprintf#59{fprintf(stderr, "%0.2lf ", R[i][j])}
      ****************************************/
      for(j = 0; j < n; j++) {
         if((i * n + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", R[i][j]);
      }
   }
   fprintf(stderr, "\nend   dump: %s\n", "R");
   fprintf(stderr, "begin dump: %s", "Q");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#71{fprintf(stderr, "\n")}
   fprintf#73{fprintf(stderr, "%0.2lf ", Q[i][j])}
   ****************************************/
   for(i = 0; i < m; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#71{fprintf(stderr, "\n")}
      fprintf#73{fprintf(stderr, "%0.2lf ", Q[i][j])}
      ****************************************/
      for(j = 0; j < n; j++) {
         if((i * n + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", Q[i][j]);
      }
   }
   fprintf(stderr, "\nend   dump: %s\n", "Q");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
/*QR Decomposition with Modified Gram Schmidt:
http://www.inf.ethz.ch/personal/gander/*/
static void kernel_gramschmidt(int m, int n, double A[1000][1200], double R[1200][1200], double Q[1000][1200]) {
   int i, j, k;
   double nrm;
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess A	 use : RW
   ****************************************/
   for(k = 0; k < n; k++) {
      nrm = 0.0;
      #pragma omp parallel for default(shared) private(i) firstprivate(m, k) reduction(+ : nrm)
      for(i = 0; i < m; i++) nrm += A[i][k] * A[i][k];
      R[k][k] = sqrt(nrm);
      #pragma omp parallel for default(shared) private(i) firstprivate(m, k)
      for(i = 0; i < m; i++) Q[i][k] = A[i][k] / R[k][k];
      #pragma omp parallel for default(shared) private(j, i) firstprivate(k, n, m)
      for(j = k + 1; j < n; j++) {
         R[k][j] = 0.0;
         // #pragma omp parallel for default(shared) private(i) firstprivate(m, k, j) reduction(+ : R[k][j])
         for(i = 0; i < m; i++) R[k][j] += Q[i][k] * A[i][j];
         // #pragma omp parallel for default(shared) private(i) firstprivate(m, k, j)
         for(i = 0; i < m; i++) A[i][j] = A[i][j] - Q[i][k] * R[k][j];
      }
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int m = 1000;
   int n = 1200;
   double (*A)[1000][1200];
   /*Variable declaration/allocation.*/
   A = (double (*)[1000][1200]) polybench_alloc_data((1000 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*R)[1200][1200];
   R = (double (*)[1200][1200]) polybench_alloc_data((1200 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*Q)[1000][1200];
   Q = (double (*)[1000][1200]) polybench_alloc_data((1000 + 0) * (1200 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(m, n, *A, *R, *Q);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_gramschmidt(m, n, *A, *R, *Q);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(m, n, *A, *R, *Q);
   /*Be clean.*/
   free((void *) A);
   ;
   free((void *) R);
   ;
   free((void *) Q);
   ;
   
   return 0;
}

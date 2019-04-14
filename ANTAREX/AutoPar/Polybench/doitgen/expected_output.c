#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "doitgen.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*doitgen.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int nr, int nq, int np, double A[150][140][160], double C4[160][160]) {
   int i, j, k;
   #pragma omp parallel for default(shared) private(i, j, k) firstprivate(nr, nq, np)
   for(i = 0; i < nr; i++) {
      // #pragma omp parallel for default(shared) private(j, k) firstprivate(nq, np, i)
      for(j = 0; j < nq; j++) {
         // #pragma omp parallel for default(shared) private(k) firstprivate(np, i, j)
         for(k = 0; k < np; k++) A[i][j][k] = (double) ((i * j + k) % np) / np;
      }
   }
   #pragma omp parallel for default(shared) private(i, j) firstprivate(np)
   for(i = 0; i < np; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(np, i)
      for(j = 0; j < np; j++) C4[i][j] = (double) (i * j % np) / np;
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int nr, int nq, int np, double A[150][140][160]) {
   int i, j, k;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "A");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#63{fprintf(stderr, "\n")}
   fprintf#65{fprintf(stderr, "%0.2lf ", A[i][j][k])}
   ****************************************/
   for(i = 0; i < nr; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#63{fprintf(stderr, "\n")}
      fprintf#65{fprintf(stderr, "%0.2lf ", A[i][j][k])}
      ****************************************/
      for(j = 0; j < nq; j++) {
         /*************** Clava msgError **************
         Variables Access as passed arguments Can not be traced inside of function calls :
         fprintf#63{fprintf(stderr, "\n")}
         fprintf#65{fprintf(stderr, "%0.2lf ", A[i][j][k])}
         ****************************************/
         for(k = 0; k < np; k++) {
            if((i * nq * np + j * np + k) % 20 == 0) fprintf(stderr, "\n");
            fprintf(stderr, "%0.2lf ", A[i][j][k]);
         }
      }
   }
   fprintf(stderr, "\nend   dump: %s\n", "A");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
void kernel_doitgen(int nr, int nq, int np, double A[150][140][160], double C4[160][160], double sum[160]) {
   int r, q, p, s;
   printf("_PB_NR = %d \t _PB_NQ = %d \t _PB_NP = %d \n", nr, nq, np);
   #pragma omp parallel for default(shared) private(r, q, p, s) firstprivate(nr, nq, np) reduction(+ : sum[:160])
   for(r = 0; r < nr; r++) {
      // #pragma omp parallel for default(shared) private(q, p, s) firstprivate(nq, np, r) reduction(+ : sum[:160])
      for(q = 0; q < nq; q++) {
         // #pragma omp parallel for default(shared) private(p, s) firstprivate(np, r, q)
         for(p = 0; p < np; p++) {
            sum[p] = 0.0;
            // #pragma omp parallel for default(shared) private(s) firstprivate(np, r, q, p) reduction(+ : sum[p])
            for(s = 0; s < np; s++) sum[p] += A[r][q][s] * C4[s][p];
         }
         // #pragma omp parallel for default(shared) private(p) firstprivate(np, r, q)
         for(p = 0; p < np; p++) A[r][q][p] = sum[p];
      }
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int nr = 150;
   int nq = 140;
   int np = 160;
   /*Variable declaration/allocation.*/   
   double (*A)[150][140][160];
   A = (double (*)[150][140][160]) polybench_alloc_data((150 + 0) * (140 + 0) * (160 + 0), sizeof(double));
   ;
   double (*sum)[160];
   sum = (double (*)[160]) polybench_alloc_data(160 + 0, sizeof(double));
   ;
   double (*C4)[160][160];
   C4 = (double (*)[160][160]) polybench_alloc_data((160 + 0) * (160 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(nr, nq, np, *A, *C4);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_doitgen(nr, nq, np, *A, *C4, *sum);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(nr, nq, np, *A);
   /*Be clean.*/
   free((void *) A);
   ;
   free((void *) sum);
   ;
   free((void *) C4);
   ;
   
   return 0;
}

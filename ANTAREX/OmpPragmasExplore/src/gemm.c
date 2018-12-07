#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "gemm.h"
#include <time.h>
#include <windows.h>
//#include <polybench.h>
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
void * polybench_alloc_data(unsigned long long n, int elt_size) {
   /// FIXME: detect overflow!
   size_t val = n;
   val *= elt_size;
   void * ret = malloc(val);
   
   return ret;
}

/*gemm.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int ni, int nj, int nk, double * alpha, double * beta, double C[1000][1100], double A[1000][1200], double B[1200][1100]) {
   int i, j;
   *alpha = 1.5;
   *beta = 1.2;
   //#pragma omp parallel for default(shared) private(i, j) firstprivate(ni, nj)
   for(i = 0; i < ni; i++) {
      //#pragma omp parallel for default(shared) private(j) firstprivate(nj, i, ni)
      for(j = 0; j < nj; j++) C[i][j] = (double) ((i * j + 1) % ni) / ni;
   }
   //#pragma omp parallel for default(shared) private(i, j) firstprivate(ni, nk)
   for(i = 0; i < ni; i++) {
      //#pragma omp parallel for default(shared) private(j) firstprivate(nk, i)
      for(j = 0; j < nk; j++) A[i][j] = (double) (i * (j + 1) % nk) / nk;
   }
   //#pragma omp parallel for default(shared) private(i, j) firstprivate(nk, nj)
   for(i = 0; i < nk; i++) {
      //#pragma omp parallel for default(shared) private(j) firstprivate(nj, i)
      for(j = 0; j < nj; j++) B[i][j] = (double) (i * (j + 2) % nj) / nj;
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int ni, int nj, double C[1000][1100]) {
   int i, j;
   fprintf((__acrt_iob_func(2)), "==BEGIN DUMP_ARRAYS==\n");
   fprintf((__acrt_iob_func(2)), "begin dump: %s", "C");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#67{fprintf(stderr, "\n")}
   fprintf#69{fprintf(stderr, "%0.2lf ", C[i][j])}
   ****************************************/
   for(i = 0; i < ni; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#67{fprintf(stderr, "\n")}
      fprintf#69{fprintf(stderr, "%0.2lf ", C[i][j])}
      ****************************************/
      for(j = 0; j < nj; j++) {
         if((i * ni + j) % 20 == 0) fprintf((__acrt_iob_func(2)), "\n");
         fprintf((__acrt_iob_func(2)), "%0.2lf ", C[i][j]);
      }
   }
   fprintf((__acrt_iob_func(2)), "\nend   dump: %s\n", "C");
   fprintf((__acrt_iob_func(2)), "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_gemm(int ni, int nj, int nk, double alpha, double beta, double C[1000][1100], double A[1000][1200], double B[1200][1100]) {
   int i, j, k;
   #pragma omp parallel for default(shared) private(i, j, k) firstprivate(ni, nj, beta, nk, alpha)
   for(i = 0; i < ni; i++) {
      #pragma omp parallel for default(shared) private(j) firstprivate(nj, i, beta)
      for(j = 0; j < nj; j++) C[i][j] *= beta;
      #pragma omp parallel for default(shared) private(k, j) firstprivate(nk, nj, alpha, i)
      for(k = 0; k < nk; k++) {
         #pragma omp parallel for default(shared) private(j) firstprivate(nj, alpha, i, k)
         for(j = 0; j < nj; j++) C[i][j] += alpha * A[i][k] * B[k][j];
      }
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int ni = 1000;
   int nj = 1100;
   int nk = 1200;
   /*Variable declaration/allocation.*/
   double alpha;
   double beta;
   double (*C)[1000][1100];
   C = (double (*)[1000][1100]) polybench_alloc_data((1000 + 0) * (1100 + 0), sizeof(double));
   ;
   double (*A)[1000][1200];
   A = (double (*)[1000][1200]) polybench_alloc_data((1000 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*B)[1200][1100];
   B = (double (*)[1200][1100]) polybench_alloc_data((1200 + 0) * (1100 + 0), sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(ni, nj, nk, &alpha, &beta, *C, *A, *B);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   printf("Started gemm\n");
   LARGE_INTEGER clava_timing_start_0, clava_timing_end_0, clava_timing_frequency_0;
   QueryPerformanceFrequency(&clava_timing_frequency_0);
   QueryPerformanceCounter(&clava_timing_start_0);
   kernel_gemm(ni, nj, nk, alpha, beta, *C, *A, *B);
   QueryPerformanceCounter(&clava_timing_end_0);
   double clava_timing_duration_0 = ((clava_timing_end_0.QuadPart-clava_timing_start_0.QuadPart) / (double)clava_timing_frequency_0.QuadPart) * (1000);
   printf("Exec time (ms):%f\n", clava_timing_duration_0);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(ni, nj, *C);
   free((void *) C);
   /*Be clean.*/
   ;
   free((void *) A);
   ;
   free((void *) B);
   ;
   
   return 0;
}

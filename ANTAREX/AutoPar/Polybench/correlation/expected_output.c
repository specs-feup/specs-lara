#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "correlation.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*correlation.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int m, int n, double * float_n, double data[1400][1200]) {
   int i, j;
   *float_n = (double) 1400;
   #pragma omp parallel for default(shared) private(i, j)
   for(i = 0; i < 1400; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(i)
      for(j = 0; j < 1200; j++) data[i][j] = (double) (i * j) / 1200 + i;
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int m, double corr[1200][1200]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "corr");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#48{fprintf(stderr, "\n")}
   fprintf#50{fprintf(stderr, "%0.2lf ", corr[i][j])}
   ****************************************/
   for(i = 0; i < m; i++) {
      /*************** Clava msgError **************
      Variables Access as passed arguments Can not be traced inside of function calls :
      fprintf#48{fprintf(stderr, "\n")}
      fprintf#50{fprintf(stderr, "%0.2lf ", corr[i][j])}
      ****************************************/
      for(j = 0; j < m; j++) {
         if((i * m + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", corr[i][j]);
      }
   }
   fprintf(stderr, "\nend   dump: %s\n", "corr");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_correlation(int m, int n, double float_n, double data[1400][1200], double corr[1200][1200], double mean[1200], double stddev[1200]) {
   int i, j, k;
   double eps = 0.1;
   #pragma omp parallel for default(shared) private(j, i) firstprivate(m, n, float_n)
   for(j = 0; j < m; j++) {
      mean[j] = 0.0;
      // #pragma omp parallel for default(shared) private(i) firstprivate(n, j) reduction(+ : mean[j])
      for(i = 0; i < n; i++) mean[j] += data[i][j];
      mean[j] /= float_n;
   }
   #pragma omp parallel for default(shared) private(j, i) firstprivate(m, n, float_n, eps)
   for(j = 0; j < m; j++) {
      stddev[j] = 0.0;
      // #pragma omp parallel for default(shared) private(i) firstprivate(n, j) reduction(+ : stddev[j])
      for(i = 0; i < n; i++) stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
      stddev[j] /= float_n;
      stddev[j] = sqrt(stddev[j]);
      stddev[j] = stddev[j] <= eps ? 1.0 : stddev[j];
   }
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n, m, float_n)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(m, i, float_n)
      for(j = 0; j < m; j++) {
         data[i][j] -= mean[j];
         data[i][j] /= sqrt(float_n) * stddev[j];
      }
   }
   #pragma omp parallel for default(shared) private(i, j, k) firstprivate(m, n)
   for(i = 0; i < m - 1; i++) {
      corr[i][i] = 1.0;
      // #pragma omp parallel for default(shared) private(j, k) firstprivate(i, m, n)
      for(j = i + 1; j < m; j++) {
         corr[i][j] = 0.0;
         // #pragma omp parallel for default(shared) private(k) firstprivate(n, i, j) reduction(+ : corr[i][j])
         for(k = 0; k < n; k++) corr[i][j] += (data[k][i] * data[k][j]);
         corr[j][i] = corr[i][j];
      }
   }
   corr[m - 1][m - 1] = 1.0;
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int n = 1400;
   int m = 1200;
   /*Variable declaration/allocation.*/
   double float_n;
   double (*data)[1400][1200];
   data = (double (*)[1400][1200]) polybench_alloc_data((1400 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*corr)[1200][1200];
   corr = (double (*)[1200][1200]) polybench_alloc_data((1200 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*mean)[1200];
   mean = (double (*)[1200]) polybench_alloc_data(1200 + 0, sizeof(double));
   ;
   double (*stddev)[1200];
   stddev = (double (*)[1200]) polybench_alloc_data(1200 + 0, sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(m, n, &float_n, *data);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_correlation(m, n, float_n, *data, *corr, *mean, *stddev);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(m, *corr);
   free((void *) data);
   /*Be clean.*/
   ;
   free((void *) corr);
   ;
   free((void *) mean);
   ;
   free((void *) stddev);
   ;
   
   return 0;
}

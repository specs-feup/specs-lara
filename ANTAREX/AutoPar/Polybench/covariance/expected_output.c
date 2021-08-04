#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "covariance.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*covariance.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int m, int n, double *float_n, double data[1400][1200]) {
   int i, j;
   *float_n = (double) n;
   for(i = 0; i < 1400; i++)
      for(j = 0; j < 1200; j++)
         data[i][j] = ((double) i * j) / 1200;
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int m, double cov[1200][1200]) {
   int i, j;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "cov");
   for(i = 0; i < m; i++)
      for(j = 0; j < m; j++) {
         if((i * m + j) % 20 == 0) fprintf(stderr, "\n");
         fprintf(stderr, "%0.2lf ", cov[i][j]);
      }
   fprintf(stderr, "\nend   dump: %s\n", "cov");
   fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

/*Main computational kernel. The whole function will be timed,
including the call and return.*/
static void kernel_covariance(int m, int n, double float_n, double data[1400][1200], double cov[1200][1200], double mean[1200]) {
   int i, j, k;
   #pragma omp parallel for default(shared) private(j, i) firstprivate(m, n, float_n, data)
   for(j = 0; j < m; j++) {
      mean[j] = 0.0;
      // #pragma omp parallel for default(shared) private(i) firstprivate(n, j, data) reduction(+ : mean[j])
      for(i = 0; i < n; i++)
         mean[j] += data[i][j];
      mean[j] /= float_n;
   }
   #pragma omp parallel for default(shared) private(i, j) firstprivate(n, m, mean)
   for(i = 0; i < n; i++) {
      // #pragma omp parallel for default(shared) private(j) firstprivate(m, i, mean)
      for(j = 0; j < m; j++)
         data[i][j] -= mean[j];
   }
   #pragma omp parallel for default(shared) private(i, j, k) firstprivate(m, n, float_n, data)
   for(i = 0; i < m; i++) {
      // #pragma omp parallel for default(shared) private(j, k) firstprivate(i, m, n, float_n, data)
      for(j = i; j < m; j++) {
         cov[i][j] = 0.0;
         // #pragma omp parallel for default(shared) private(k) firstprivate(n, i, j, data) reduction(+ : cov[i][j])
         for(k = 0; k < n; k++)
            cov[i][j] += data[k][i] * data[k][j];
         cov[i][j] /= (float_n - 1.0);
         cov[j][i] = cov[i][j];
      }
   }
}

int main(int argc, char **argv) {
   /*Retrieve problem size.*/
   int n = 1400;
   int m = 1200;
   /*Variable declaration/allocation.*/
   double float_n;
   double (*data)[1400][1200];
   data = (double (*)[1400][1200]) polybench_alloc_data((1400 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*cov)[1200][1200];
   cov = (double (*)[1200][1200]) polybench_alloc_data((1200 + 0) * (1200 + 0), sizeof(double));
   ;
   double (*mean)[1200];
   mean = (double (*)[1200]) polybench_alloc_data(1200 + 0, sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(m, n, &float_n, *data);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_covariance(m, n, float_n, *data, *cov, *mean);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(m, *cov);
   /*Be clean.*/
   free((void *) data);
   ;
   free((void *) cov);
   ;
   free((void *) mean);
   ;
   
   return 0;
}

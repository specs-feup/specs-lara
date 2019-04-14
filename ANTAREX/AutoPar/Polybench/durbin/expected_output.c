#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <polybench.h>
#include "durbin.h"
/**
* This version is stamped on May 10, 2016
*
* Contact:
*   Louis-Noel Pouchet <pouchet.ohio-state.edu>
*   Tomofumi Yuki <tomofumi.yuki.fr>
*
* Web address: http://polybench.sourceforge.net
*/
/*durbin.c: this file is part of PolyBench/C*/
/*Include polybench common header.*/
/*Include benchmark-specific header.*/
/*Array initialization.*/
static void init_array(int n, double r[2000]) {
   int i, j;
   #pragma omp parallel for default(shared) private(i) firstprivate(n)
   for(i = 0; i < n; i++) {
      r[i] = (n + 1 - i);
   }
}

/*DCE code. Must scan the entire live-out data.
Can be used also to check the correctness of the output.*/
static void print_array(int n, double y[2000]) {
   int i;
   fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
   fprintf(stderr, "begin dump: %s", "y");
   /*************** Clava msgError **************
   Variables Access as passed arguments Can not be traced inside of function calls :
   fprintf#40{fprintf(stderr, "\n")}
   fprintf#42{fprintf(stderr, "%0.2lf ", y[i])}
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
static void kernel_durbin(int n, double r[2000], double y[2000]) {
   double z[2000];
   double alpha;
   double beta;
   double sum;
   int i, k;
   y[0] = -r[0];
   beta = 1.0;
   alpha = -r[0];
   /*************** Clava msgError **************
   Variable alpha could not be categorized into any OpenMP Variable Scopeuse : RWR
   Variable beta could not be categorized into any OpenMP Variable Scopeuse : RWR
   ****************************************/
   for(k = 1; k < n; k++) {
      beta = (1 - alpha * alpha) * beta;
      sum = 0.0;
      #pragma omp parallel for default(shared) private(i) firstprivate(k) reduction(+ : sum)
      for(i = 0; i < k; i++) {
         sum += r[k - i - 1] * y[i];
      }
      alpha = -(r[k] + sum) / beta;
      #pragma omp parallel for default(shared) private(i) firstprivate(k, alpha)
      for(i = 0; i < k; i++) {
         z[i] = y[i] + alpha * y[k - i - 1];
      }
      #pragma omp parallel for default(shared) private(i) firstprivate(k)
      for(i = 0; i < k; i++) {
         y[i] = z[i];
      }
      y[k] = alpha;
   }
}

int main(int argc, char ** argv) {
   /*Retrieve problem size.*/
   int n = 2000;
   double (*r)[2000];
   /*Variable declaration/allocation.*/
   r = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   double (*y)[2000];
   y = (double (*)[2000]) polybench_alloc_data(2000 + 0, sizeof(double));
   ;
   /*Initialize array(s).*/
   init_array(n, *r);
   /*Start timer.*/
   ;
   /*Run kernel.*/
   kernel_durbin(n, *r, *y);
   /*Stop and print timer.*/
   ;
   ;
   /*Prevent dead-code elimination. All live-out data must be printed
   by the function call in argument.*/
   if(argc > 42 && !strcmp(argv[0], "")) print_array(n, *y);
   /*Be clean.*/
   free((void *) r);
   ;
   free((void *) y);
   ;
   
   return 0;
}

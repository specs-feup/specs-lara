#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
/*
Simple matrix multiplication example.
*/
/*
matrix multiplication
*/

void matrix_mult(int const N, int const M, int const K, double A[N][M], double B[M][K], double C[N][K]) {
   for(int i = 0; i < N; i++) {
      for(int j = 0; j < K; j++) {
         C[i][j] = 0;
      }
   }
   #pragma loop1
   #pragma omp parallel for default(shared) firstprivate(N, K, M, A, B)
   for(int i = 0; i < N; i++) {
      for(int j = 0; j < K; j++) {
         for(int l = 0; l < M; l++) {
            C[i][j] += A[i][l] * B[l][j];
         }
      }
   }
}

/*
* Set an N by M matrix A to random values
*/

void init_matrix(int const N, int const M, double A[N][M]) {
   for(int i = 0; i < N; ++i) {
      for(int j = 0; j < M; ++j) {
         A[i][j] = ((double) rand()) / (double) 2147483647;
      }
   }
}


void print_matrix_result(int const N, int const K, double A[N][K]) {
   double acc = 0.0;
   for(int i = 0; i < N; ++i) {
      for(int j = 0; j < K; ++j) {
         acc += A[i][j];
      }
   }
   printf("Result acc: %f\n", acc);
}


void test_matrix_mul() {
   /*
   int N=2048;
   int M=1024;
   int K=2048;
   */
   int N = 1024;
   int M = 1024;
   int K = 512;
   // allocate matrices
   double (*A)[M] = malloc(sizeof(double[N][M]));
   double (*B)[K] = malloc(sizeof(double[M][K]));
   double (*C)[K] = malloc(sizeof(double[N][K]));
   // initialize matrices
   init_matrix(N, M, A);
   init_matrix(M, K, B);
   // do: C = A*B
   matrix_mult(N, M, K, A, B, C);
   print_matrix_result(N, K, C);
   // free	
   free(C);
   free(B);
   free(A);
}


int main() {
   // To make results repeatable
   srand(0);
   test_matrix_mul();
}

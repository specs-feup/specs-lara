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

void matrix_mult(double const * A, double const * B, double * C, int const N, int const M, int const K)
{
   int BS2 = 64;
   int BS1 = 64;
   for (int ii = 0; ii < N; ii++)
   {
      for (int jj = 0; jj < K; jj++)
      {
         //C[i][j] = 0;
         C[K * ii + jj] = 0;
      }
   }
   for (size_t l_block = 0; l_block < M; l_block = l_block + (BS1))
   {
      for (size_t j_block = 0; j_block < K; j_block = j_block + (BS2))
      {
         for (int i = 0; i < N; i++)
         {
            int l_limit = l_block + BS1;
            if (l_limit > M) l_limit = M;
            for (int l = l_block; l < l_limit; l++)
            {
               int j_limit = j_block + BS2;
               if (j_limit > K) j_limit = K;
               for (int j = j_block; j < j_limit; j++)
               {
                  //C[i][j] += A[i][l]*B[l][j];
                  C[K * i + j] = C[K * i + j] + (A[M * i + l] * B[K * l + j]);
               }
            }
         }
      }
   }
}

void matrix_mult_call_specialization(double A[131072], double B[131072], double C[262144], int const N, int const M, int const K)
{
   int BS2 = 64;
   int BS1 = 64;
   for (int ii = 0; ii < N; ii++)
   {
      for (int jj = 0; jj < K; jj++)
      {
         C[K * ii + jj] = 0;
      }
   }
   for (int l_block = 0; l_block < M; l_block = l_block + (BS1))
   {
      for (int j_block = 0; j_block < K; j_block = j_block + (BS2))
      {

         for (int i = 0; i < N; i++)
         {
            int l_limit = l_block + BS1;
            if (l_limit > M) l_limit = M;
            for (int l = l_block; l < l_limit; l++)
            {
               int j_limit = j_block + BS2;
               if (j_limit > K) j_limit = K;
               for (int j = j_block; j < j_limit; j++)
               {
                  C[K * i + j] = C[K * i + j] + (A[M * i + l] * B[K * l + j]);
               }
            }
         }
      }
   }
}

void matrix_mult_call_specialization_0(double A[131072], double B[131072], double C[262144], int const N, int const M, int const K)
{
   int BS2 = 64;
   int BS1 = 64;
   /*************** Clava msgError **************
         output   C[K * ii + jj]#45 -> C[K * ii + jj]#45 +,*   MoZ   sameloop IsdependentOuterloop=true  IsdependentCurrentloop=true   IsdependentInnerloop=true
         check file : /deploop#42[matrix_mult_call_specialization_0]Array_C.t
   ****************************************/
   for (int ii = 0; ii < N; ii++)
   {
      #pragma omp parallel for default(shared) firstprivate(C, K, ii)
      // deploop#43[matrix_mult_call_specialization_0]Array_C.t
      for (int jj = 0; jj < K; jj++)
      {
         //C[i][j] = 0;
         C[K * ii + jj] = 0;
      }
   }
   #pragma omp parallel for default(shared) firstprivate(A, B, M, BS1, K, BS2, N) reduction (+:C[:262144])
   // deploop#48[matrix_mult_call_specialization_0]Array_A_B_C.t
   for (int l_block = 0; l_block < M; l_block = l_block + (BS1))
   {
      //#pragma omp parallel for default(shared) firstprivate(A, B, K, BS2, N, l_block, BS1, M) reduction (+:C[:262144])
      // deploop#49[matrix_mult_call_specialization_0]Array_A_B_C.t
      for (int j_block = 0; j_block < K; j_block = j_block + (BS2))
      {
         /*************** Clava msgError **************
               anti  C[K * i + j]#58 -> C[K * i + j]#58  0,+,* MoZ   sameloop IsdependentOuterloop=true  IsdependentCurrentloop=true   IsdependentInnerloop=true
               anti  C[K * i + j]#58 -> C[K * i + j]#58  +,*,* MoZ   sameloop IsdependentOuterloop=true  IsdependentCurrentloop=true   IsdependentInnerloop=true
               flow  C[K * i + j]#58 -> C[K * i + j]#58  0,+,* MoZ   sameloop IsdependentOuterloop=true  IsdependentCurrentloop=true   IsdependentInnerloop=true
               flow  C[K * i + j]#58 -> C[K * i + j]#58  +,*,* MoZ   sameloop IsdependentOuterloop=true  IsdependentCurrentloop=true   IsdependentInnerloop=true
               output   C[K * i + j]#58 -> C[K * i + j]#58  +,*,* MoZ   sameloop IsdependentOuterloop=true  IsdependentCurrentloop=true   IsdependentInnerloop=true
               check file : /deploop#50[matrix_mult_call_specialization_0]Array_A_B_C.t
         ****************************************/
         for (int i = 0; i < N; i++)
         {
            int l_limit = l_block + BS1;
            if (l_limit > M) l_limit = M;
            //#pragma omp parallel for default(shared) firstprivate(A, B, l_block, l_limit, j_block, BS2, K, M, i) reduction (+:C[:262144])
            // deploop#53[matrix_mult_call_specialization_0]Array_A_B_C.t
            for (int l = l_block; l < l_limit; l++)
            {
               int j_limit = j_block + BS2;
               if (j_limit > K) j_limit = K;
               //#pragma omp parallel for default(shared) firstprivate(A, B, C, j_block, j_limit, M, i, K, l)
               // deploop#56[matrix_mult_call_specialization_0]Array_A_B_C.t
               for (int j = j_block; j < j_limit; j++)
               {
                  //C[i][j] += A[i][l]*B[l][j];
                  C[K * i + j] = C[K * i + j] + (A[M * i + l] * B[K * l + j]);
               }
            }
         }
      }
   }
}

/*
* Set an N by M matrix A to random values
*/

void init_matrix(double * A, int const N, int const M)
{
   for (int i = 0; i < N; ++i)
   {
      for (int j = 0; j < M; ++j)
      {
         //A[i][j] = ((double) rand()) / (double) RAND_MAX;
         A[M * i + j] = ((double) rand()) / (double) 32767;
      }
   }
}


void print_matrix_result(double * A, int const N, int const K)
{
   double acc = 0.0;
   for (int i = 0; i < N; ++i)
   {
      for (int j = 0; j < K; ++j)
      {
         //acc += A[i][j];
         acc = acc + (A[K * i + j]);
      }
   }
   printf("Result acc: %f\n", acc);
}


void test_matrix_mul()
{
   int N = 512;
   int M = 256;
   int K = 512;
   //double A[N][M];
   //double B[M][K];
   //double C[N][K];
   double * A = (double *) malloc(N * M * sizeof(double));
   double * B = (double *) malloc(M * K * sizeof(double));
   double * C = (double *) malloc(N * K * sizeof(double));
   double * C_OMP = (double *) malloc(N * K * sizeof(double));
   // initialize matrices
   init_matrix(A, N, M);
   init_matrix(B, M, K);
   // do: C = A*B
   matrix_mult_call_specialization(A, B, C, N, M, K);
   matrix_mult_call_specialization_0(A, B, C_OMP, N, M, K);
   print_matrix_result(C, N, K);
   print_matrix_result(C_OMP, N, K);
   free(A);
   free(B);
   free(C);
}


int main()
{
   // To make results repeatable
   srand(0);
   test_matrix_mul();
}


/*
	Simple matrix multiplication example.
*/

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>

/*
	matrix multiplication
*/
void matrix_mult(const double* A , const double* B, double* C, const int N, const int M, const int K) {

   for(int ii=0; ii<N; ii++) {
       for(int jj=0; jj<K; jj++) {
           C[K*ii + jj] = 0;
       }
    }

    // main loop nest
    for(int i=0; i<N; i++) {
        for(int l=0; l< M; l++) {
            for(int j=0; j< K; j++) {
                C[K*i + j] += A[M*i+l]*B[K*l+j];
            }
        }
    }
}

/*
 * Set an N by M matrix A to random values
 */
void init_matrix(double *A, const int N, const int M) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
	       //A[i][j] = ((double) rand()) / (double) RAND_MAX; 
	       A[M*i + j] = ((double) rand()) / (double) RAND_MAX; 			
		}
	}
}

void print_matrix_result(double *A, const int N, const int K) {
	double acc = 0.0;
	
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < K; ++j) {
			//acc += A[i][j];
			acc += A[K*i + j];
		}
	}
		
		

	printf("Result acc: %f\n", acc);
}

	
	
void test_matrix_mul() {	
	int N=512; 
	int M=256;
	int K=512;
	
	//double A[N][M];
	//double B[M][K];
	//double C[N][K];
	double* A = (double *) malloc(N*M*sizeof(double)); 
	double* B = (double *) malloc(M*K*sizeof(double));
	double* C = (double *) malloc(N*K*sizeof(double));
		
	// initialize matrices
	init_matrix(A, N, M);
	init_matrix(B, M, K);
	
	#pragma clava opencl_call matrix_opencl.json
	matrix_mult(A, B, C, N, M, K);
	
	print_matrix_result(C, N, K);
	
}

int main() {
	
	// To make results repeatable
	srand(0);

	test_matrix_mul();	
}	

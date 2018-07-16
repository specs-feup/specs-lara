/*
	Simple matrix multiplication example.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
	matrix multiplication
*/
void matrix_mult(const int N, const int M, const int K, double A[N][M], double B[M][K], double C[N][K]) {

   for(int i=0; i<N; i++) {
       for(int j=0; j<K; j++) {
           C[i][j] = 0;
       }
    }

	#pragma loop1
    for(int i=0; i<N; i++) {
        for(int l=0; l< M; l++) {
            for(int j=0; j< K; j++) {
				C[i][j] += A[i][l]*B[l][j];
            }
        }
    }
}

/*
 * Set an N by M matrix A to random values
 */
void init_matrix(const int N, const int M, double A[N][M]) {
	
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < M; ++j) {
			
			A[i][j] = ((double) rand()) / (double) RAND_MAX;	
		}
	}
}

void print_matrix_result(const int N, const int K, double A[N][K]) {
	
	double acc = 0.0;
	
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < K; ++j) {
			
			acc += A[i][j];
		}
	}
		
	printf("Result acc: %f\n", acc);
}

double** alloc_matrix(const int N, const int M) {
	
	double* AI = malloc(N*M*sizeof(*AI));
	double** A = malloc(N*sizeof(*A));
	
	for(int i = 0; i < N; i++) {
		A[i] = &(AI[i * M]);
	}
	
	return A;
}

void free_matrix(double** A) {
	
	free(A[0]);
	free(A);
}
	
	
void test_matrix_mul() {	
	
	int N=2048; 
	int M=1024;
	int K=2048;
	
	double (*A)[M] = malloc(sizeof(double[N][M]));
	double (*B)[K] = malloc(sizeof(double[M][K]));
	double (*C)[K] = malloc(sizeof(double[N][K]));
	
	//~ double** A = alloc_matrix(N, M);
	//~ double** B = alloc_matrix(M, K);
	//~ double** C = alloc_matrix(N, K);
			
	// initialize matrices
	init_matrix(N, M, A);
	init_matrix(M, K, B);
	
	// do: C = A*B
	matrix_mult(N, M, K, A, B, C);
	
	print_matrix_result(N, K, C);
		
	// free
	//~ free_matrix(C);
	//~ free_matrix(B);
	//~ free_matrix(A);
	
	free(C);
	free(B);
	free(A);
}

int main() {
	
	// To make results repeatable
	srand(0);

	test_matrix_mul();	
}	

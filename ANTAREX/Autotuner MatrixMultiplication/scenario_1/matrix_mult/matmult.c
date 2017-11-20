/*
	Simple matrix multiplication example.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <margot.h>

#define SCENARIO 1 // 1: one execution; 2: N_EXEC executions with different matrix sizes

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef double REAL;

/*
	matrix multiplication with loop tiling
*/
void matrix_mult_tiling(const REAL* A , const REAL* B, REAL* C, const int N, const int M, const int K, const int BS1, const int BS2) {

   for(int i=0; i<N; i++) {
       for(int j=0; j<K; j++) {
           C[K*i + j] = 0;
       }
    }

    for(int l2=0; l2<M; l2 += BS1) {
        for(int j2=0; j2<K; j2 += BS2) {
            for(int i=0; i<N; i++) {
                for(int l=l2; l< MIN(M, l2+BS1); l++) {
                    for(int j=j2; j< MIN(K, j2+BS2); j++) {
                        C[K*i + j] += A[M*i+l]*B[K*l+j];
                    }
                }
            }
        }
    }
}

/*
 * Set an N by M matrix A to random values
 */
void init_matrix(REAL *A, const int N, const int M) {
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j)
	       A[M*i + j] = (REAL) ((double) rand()) / (double) RAND_MAX;
}

/*
	Two possible execution scenarios:
*/

#if SCENARIO == 1

int main() {

  // margot initialize the framework
  margot_init();

	int BS1 = 64; // tiling sizes: 8, 16, 32, 64
	int BS2 = 64; // tiling size: 8, 16, 32, 64

	int N=512;
	int M=256;
	int K=512;

	REAL* A = (REAL *) malloc(N*M*sizeof(REAL));
	REAL* B = (REAL *) malloc(M*K*sizeof(REAL));
	REAL* C = (REAL *) malloc(N*K*sizeof(REAL));

	// initialize matrices
	init_matrix(A, N, M);
	init_matrix(B, M, K);

  // margot wrapper code begin
  if (margot_matmul_update( &BS1, &BS2 ))
  {
    margot_matmul_configuration_applied();
  }
  margot_matmul_start_monitor();

	// do: C = A*B
	matrix_mult_tiling(A , B, C, N, M, K, BS1, BS2);

  // margot wrapper code end
  margot_matmul_stop_monitor();
  margot_matmul_log();

}

#else

#define N_EXEC 10
int N1[N_EXEC] = {512, 512, 512, 256, 256, 512, 128, 128, 256, 512};
int M1[N_EXEC] = {512, 512, 256, 256, 256, 512, 128, 256, 256, 512};
int K1[N_EXEC] = {512, 256, 512, 256, 128, 512, 128, 128, 256, 512};

int main() {

  // margot initialize the framework
  margot_init();

	int BS1 = 64; // tiling sizes: 8, 16, 32, 64
	int BS2 = 64; // tiling size: 8, 16, 32, 64

	int N;
	int M;
	int K;

	REAL* A;
	REAL* B;
	REAL* C;

	for (int i=0; i<N_EXEC; i++) {

		N= N1[i];
		M= M1[i];
		K= K1[i];

		A = (REAL *) malloc(N*M*sizeof(REAL));
		B = (REAL *) malloc(M*K*sizeof(REAL));
		C = (REAL *) malloc(N*K*sizeof(REAL));

		// initialize matrices
		init_matrix(A, N, M);
		init_matrix(B, M, K);

    // margot wrapper code begin
    if (margot_matmul_update( &BS1, &BS2 ))
    {
      margot_matmul_configuration_applied();
    }
    margot_matmul_start_monitor();

		// do: C = A*B
		matrix_mult_tiling(A , B, C, N, M, K, BS1, BS2);

    // margot wrapper code end
    margot_matmul_stop_monitor();
    margot_matmul_log();
	}
}
#endif

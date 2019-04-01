#include "matrix_mul.hpp"

template< typename T >
void matrix_mult(const std::vector<T>& A , const std::vector<T>& B, std::vector<T>& C, const int N, const int M, const int K) {

   for(int i=0; i<N; i++) {
       for(int j=0; j<K; j++) {
           C[K*i + j] = 0;
       }
    }

	for(int i=0; i<N; i++) {
		for(int l=0; l < M; l++) {
			for(int j=0; j < K; j++) {
				C[K*i + j] += A[M*i+l]*B[K*l+j];
			}
		}
	}

}

int main() {
  // iteration counter for printing purpose
  int num_iters = 21;

  // matrix sizes
  int N,M,K;

  // this is the main loop of the application
  for( int i = 0; i < num_iters; i++ )
  {

    // declare the matrices
    std::vector<double> A,B,C;

    // mimic the read of the input matrices
    N = N_sizes[i];
    M = M_sizes[i];
    K = K_sizes[i];
    init_matrix(A, N, M);
    init_matrix(B, M, K);
    init_matrix(C, N, K, false);

    // execute the kernel
	std::chrono::high_resolution_clock::time_point timing_start = std::chrono::high_resolution_clock::now();
	matrix_mult(A, B, C, N, M, K);
	std::chrono::high_resolution_clock::time_point timing_end = std::chrono::high_resolution_clock::now();
	auto timing_duration = std::chrono::duration_cast<std::chrono::milliseconds>(timing_end - timing_start).count();
	std::cout << timing_duration << "ms" << std::endl;

    std::cout << "#" << i << " ";
    std::cout << "C[0][0] = " << C[0] << " ";
    std::cout << "[" << N << "x" << M << "]";
    std::cout << " X ";
    std::cout << "[" << M << "x" << K << "]" << std::endl;
  }


  return EXIT_SUCCESS;
}

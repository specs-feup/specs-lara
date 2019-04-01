#include "matrix_mul.hpp"
#include <margot.hpp>
template <typename T>
void matrix_mult_tiling(std::vector<T> const& A, std::vector<T> const& B, std::vector<T>& C, int const N, int const M, int const K, int const BS1, int const BS2) {
   
   for(int i = 0; i < N; i++) {
      for(int j = 0; j < K; j++) {
         C[K * i + j] = 0;
      }
   }
   
   for(int l2 = 0; l2 < M; l2 += BS1) {
      for(int j2 = 0; j2 < K; j2 += BS2) {
         for(int i = 0; i < N; i++) {
            for(int l = l2; l < std::min(M, l2 + BS1); l++) {
               for(int j = j2; j < std::min(K, j2 + BS2); j++) {
                  C[K * i + j] += A[M * i + l] * B[K * l + j];
               }
            }
         }
      }
   }
}

int main() {
   margot::init();
  // iteration counter for printing purpose
  int num_iters = 21;

   // default value of the parameters
   int BS1 = 64;
   int BS2 = 64;

  // matrix sizes
  int N,M,K;

  // this is the main loop of the application
  for( int i = 0; i < 1; i++ )
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
	if(margot::matmul::update(BS1,BS2, N, M, K)) {
		margot::matmul::manager.configuration_applied();
	}
	margot::matmul::start_monitor();
	matrix_mult_tiling(A, B, C, N, M, K, BS1, BS2); // original call
	margot::matmul::stop_monitor();
	margot::matmul::log();


    std::cout << "#" << i << " ";
    std::cout << "C[0][0] = " << C[0] << " ";
    std::cout << "[" << N << "x" << M << "]";
    std::cout << " X ";
    std::cout << "[" << M << "x" << K << "]" << std::endl;
  }


  return EXIT_SUCCESS;
}

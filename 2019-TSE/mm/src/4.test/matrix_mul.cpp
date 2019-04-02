#include "matrix_mul.hpp"
#include "rapl.h"

template <typename T>
void matrix_mult_tiling(std::vector<T> const& A, std::vector<T> const& B, std::vector<T>& C, int const N, int const M, int const K, int const BS1, int const BS2) {
   
   for(int i = 0; i < N; i++) {
      for(int j = 0; j < K; j++) {
         C[K * i + j] = 0;
      }
   }
   
   for(int l2 = 0; l2 < M; l2 += BS1) {
	   int l_bound = std::min(M, l2 + BS1);
      for(int j2 = 0; j2 < K; j2 += BS2) {
		  int j_bound = std::min(K, j2 + BS2);
         for(int i = 0; i < N; i++) {
            for(int l = l2; l < l_bound; l++) {
               for(int j = j2; j < j_bound; j++) {
                  C[K * i + j] += A[M * i + l] * B[K * l + j];
               }
            }
         }
      }
   }
}

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

int main(int argc, char** argv) {
	
	int runs = std::stoi(argv[1]);
	int size = std::stoi(argv[2]);
	int tile = std::stoi(argv[3]);

	printf("RUNS=%d, SIZE=%d", runs, size);
	if(tile != 0) {
		printf(", TILE=%d", tile);
	}
	printf("\n");

  // matrix sizes
  int N,M,K;

	long long total_energy = 0;
	auto total_time = 0;

  // this is the main loop of the application
  for( int i = 0; i < runs; i++ )
  {

    // declare the matrices
    std::vector<double> A,B,C;

    // mimic the read of the input matrices
    N = size;
    M = size;
    K = size;
    init_matrix(A, N, M);
    init_matrix(B, M, K);
    init_matrix(C, N, K, false);

    // execute the kernel
	if(tile != 0) {
		auto e0 = rapl_energy();
		std::chrono::high_resolution_clock::time_point timing_start = std::chrono::high_resolution_clock::now();
		matrix_mult_tiling(A, B, C, N, M, K, tile, tile);
		std::chrono::high_resolution_clock::time_point timing_end = std::chrono::high_resolution_clock::now(); 
		auto e1 = rapl_energy();
		total_energy += (e1-e0);
		total_time += std::chrono::duration_cast<std::chrono::milliseconds>(timing_end - timing_start).count();
	} else {
		auto e0 = rapl_energy();
		std::chrono::high_resolution_clock::time_point timing_start = std::chrono::high_resolution_clock::now();
		matrix_mult(A, B, C, N, M, K);
		std::chrono::high_resolution_clock::time_point timing_end = std::chrono::high_resolution_clock::now(); 
		auto e1 = rapl_energy();
		total_energy += (e1-e0);
		total_time += std::chrono::duration_cast<std::chrono::milliseconds>(timing_end - timing_start).count();
	}

    std::cout << "#" << i << " ";
    std::cout << "C[0][0] = " << C[0] << " ";
    std::cout << "[" << N << "x" << M << "]";
    std::cout << " X ";
    std::cout << "[" << M << "x" << K << "]" << std::endl;
  }
	std::cout << 1.0 * total_time / runs << " ms" << std::endl;
	std::cout << 1.0 * total_energy / runs << " uJ" << std::endl;

  return EXIT_SUCCESS;
}

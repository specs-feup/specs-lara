#include "matrix_mul.hpp"

template< typename T >
void matrix_mult_tiling(const std::vector<T>& A , const std::vector<T>& B, std::vector<T>& C, const int N, const int M, const int K, const int BS1, const int BS2) {

   for(int i=0; i<N; i++) {
       for(int j=0; j<K; j++) {
           C[K*i + j] = 0;
       }
    }

    for(int l2=0; l2<M; l2 += BS1) {
        for(int j2=0; j2<K; j2 += BS2) {
            for(int i=0; i<N; i++) {
                for(int l=l2; l< std::min(M, l2+BS1); l++) {
                    for(int j=j2; j< std::min(K, j2+BS2); j++) {
                        C[K*i + j] += A[M*i+l]*B[K*l+j];
                    }
                }
            }
        }
    }
}

int main() {
  // iteration counter for printing purpose
  std::size_t number_of_the_iteration = 0;

  // default value of the parameters
  int BS1 = 64;
  int BS2 = 64;

  // data features
  int N,M,K;

  // the object that mimics reading a matrix
  struct random_size_generator_t get_size;

  // define the time of execution
  const auto time_to_execute = std::chrono::seconds(5);
  const auto stop_time_point = std::chrono::steady_clock::now() + time_to_execute;

  // this is the main loop of the application
  while( std::chrono::steady_clock::now() < stop_time_point )
  {

    // declare the matrices
    std::vector<double> A,B,C;

    // mimic the read of the input matrices
    N = get_size();
    M = get_size();
    K = get_size();
    init_matrix(A, N, M);
    init_matrix(B, M, K);
    init_matrix(C, N, K, false);

    // execute the kernel
    matrix_mult_tiling(A, B, C, N, M, K, BS1, BS2);

    std::cout << "#" << number_of_the_iteration++ << " ";
    std::cout << "[" << N << "x" << M << "]";
    std::cout << " X ";
    std::cout << "[" << M << "x" << K << "]" << " ";
    std::cout << "   BS1: " << BS1 << " BS2: " << BS2 << std::endl;
  }


  return EXIT_SUCCESS;
}

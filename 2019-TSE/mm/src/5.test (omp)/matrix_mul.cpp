#include "matrix_mul.hpp"
#include <fstream>

template <typename T>
void matrix_mult_tiling(std::vector<T> const& A, std::vector<T> const& B, std::vector<T>& C, int const N, int const M, int const K, int const BS1, int const BS2, int const BS3)
{

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < K; j++) {
            C[K * i + j] = 0;
        }
    }

    for(int i2= 0; i2 < N; i2 += BS3) {
        int i_bound = std::min(N, i2 + BS3);
        for(int l2 = 0; l2 < M; l2 += BS1) {
            int l_bound = std::min(M, l2 + BS1);
            for(int j2 = 0; j2 < K; j2 += BS2) {
                int j_bound = std::min(K, j2 + BS2);
                for(int i = i2; i < i_bound; i++) {
                    for(int l = l2; l < l_bound; l++) {
                        for(int j = j2; j < j_bound; j++) {
                            C[K * i + j] += A[M * i + l] * B[K * l + j];
                        }
                    }
                }
            }
        }
    }
}

template< typename T >
void matrix_mult(const std::vector<T>& A, const std::vector<T>& B, std::vector<T>& C, const int N, const int M, const int K)
{

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

int main(int argc, char** argv)
{
    int size = 2048;
    int tile = 1024;

    // matrix sizes
    int N,M,K;

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
	matrix_mult_tiling(A, B, C, N, M, K, tile, tile, tile);
	
	std::ofstream out("output.txt");
	out << C[0];

    return EXIT_SUCCESS;
}

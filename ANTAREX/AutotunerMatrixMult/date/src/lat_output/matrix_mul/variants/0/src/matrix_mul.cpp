#include <chrono>
#include <fstream>
#include <iostream>
#include "matrix_mul.hpp"
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
   std::ofstream log_file_0;
   log_file_0.open("/home/pedro/Documents/repositories/specs-lara/ANTAREX/AutotunerMatrixMult/date/src/lat_output/matrix_mul/results/times.txt", std::ios_base::app);
   // iteration counter for printing purpose
   std::size_t number_of_the_iteration = 0;
   // default value of the parameters
   int BS1 = 64;
   int BS2 = 64;
   // data features
   int N;
   int M;
   int K;
   // the object that mimics reading a matrix
   struct random_size_generator_t get_size;
   // define the time of execution
   auto const time_to_execute = std::chrono::seconds(100);
   auto const stop_time_point = std::chrono::steady_clock::now() + time_to_execute;
   // this is the main loop of the application
   while(std::chrono::steady_clock::now() < stop_time_point) {
      // declare the matrices
      std::vector<double> A;
      std::vector<double> B;
      std::vector<double> C;
      // mimic the read of the input matrices
      N = get_size();
      M = get_size();
      K = get_size();
      init_matrix(A, N, M);
      init_matrix(B, M, K);
      init_matrix(C, N, K, false);
      // execute the kernel
      std::chrono::high_resolution_clock::time_point clava_timing_start_0 = std::chrono::high_resolution_clock::now();
      std::chrono::high_resolution_clock::time_point clava_timing_start_1 = std::chrono::high_resolution_clock::now();
      matrix_mult_tiling(A, B, C, N, M, K, 8, 8);
      std::chrono::high_resolution_clock::time_point clava_timing_end_1 = std::chrono::high_resolution_clock::now();
      auto clava_timing_duration_1 = std::chrono::duration_cast<std::chrono::nanoseconds>(clava_timing_end_1 - clava_timing_start_1).count();
      std::cout << "elapsed:" << clava_timing_duration_1 << std::endl;
      std::chrono::high_resolution_clock::time_point clava_timing_end_0 = std::chrono::high_resolution_clock::now();
      auto clava_timing_duration_0 = std::chrono::duration_cast<std::chrono::nanoseconds>(clava_timing_end_0 - clava_timing_start_0).count();
      log_file_0 << "matrix_mul_0 0 " << clava_timing_duration_0 << "ns" << std::endl;
      std::cout << "#" << number_of_the_iteration++ << " ";
      std::cout << "[" << N << "x" << M << "]";
      std::cout << " X ";
      std::cout << "[" << M << "x" << K << "]" << " ";
      std::cout << "   BS1: " << 8 << " BS2: " << 8 << std::endl;
   }
   
   return 0;
}


void matrix_mult_tiling(std::vector<double> const& A, std::vector<double> const& B, std::vector<double>& C, int const N, int const M, int const K, int const BS1, int const BS2) {
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

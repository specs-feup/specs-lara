#ifndef _MATRIX_MUL_HPP_
#define _MATRIX_MUL_HPP_

#include <iostream>
#include <chrono>
#include <thread>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <array>
#include <vector>
template <typename T>
void init_matrix(std::vector<T>& A, int const N, int const M, bool const initalize = true) {
   A.resize(N * M);
   if(!initalize) {
      
      return;
   }
   std::default_random_engine generator;
   std::uniform_real_distribution<T> distribution((0.0, 1.0));
   for(int i = 0; i < N; ++i) for(int j = 0; j < M; ++j) A[M * i + j] = distribution(generator);
}

//~ template< typename T >
//~ void matrix_mult_tiling(const std::vector<T>& A , const std::vector<T>& B, std::vector<T>& C, const int N, const int M, const int K, const int BS1, const int BS2);

struct random_size_generator_t {
   static int const number_of_sizes = 3;
   static std::array<int, number_of_sizes> const SPACE{{512, 256, 128}};
   std::default_random_engine generator;
   std::uniform_int_distribution<int> distribution;
   random_size_generator_t() : distribution(std::uniform_int_distribution<int>(0, number_of_sizes - 1)) {
   }
   
   inline int operator()() {
      
      return SPACE[this->distribution(this->generator)];
   }
};

std::array<int, random_size_generator_t::number_of_sizes> const random_size_generator_t::SPACE;

void init_matrix(std::vector<double>& A, int const N, int const M, bool const initalize = true) {
   A.resize(N * M);
   if(!initalize) {
      
      return;
   }
   std::default_random_engine generator;
   std::uniform_real_distribution<double> distribution(0.0, 1.0);
   for(int i = 0; i < N; ++i) for(int j = 0; j < M; ++j) A[M * i + j] = distribution(generator);
}

#endif

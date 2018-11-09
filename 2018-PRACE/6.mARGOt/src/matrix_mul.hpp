#include <iostream>
#include <chrono>
#include <thread>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <array>
#include <vector>

template<typename T>
void init_matrix(std::vector<T>& A, const int N, const int M, const bool initalize = true )
{
  A.resize(N*M);
  if (!initalize)
  {
    return;
  }
  std::default_random_engine generator;
  std::uniform_real_distribution<T> distribution(0.0,1.0);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < M; ++j)
      A[M*i + j] = distribution(generator);
}

struct random_size_generator_t
{
  static constexpr int number_of_sizes = 3;
  static constexpr std::array<int, number_of_sizes> SPACE {{512, 256, 128}};
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution;


  random_size_generator_t( void )
    : distribution(std::uniform_int_distribution<int>(0,number_of_sizes - 1))
  {}

  inline int operator()( void )
  {
    return SPACE[distribution(generator)];
  }
};

constexpr std::array<int, random_size_generator_t::number_of_sizes> random_size_generator_t::SPACE;

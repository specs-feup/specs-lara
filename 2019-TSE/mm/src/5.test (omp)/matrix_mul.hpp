#include <iostream>
#include <chrono>
#include <thread>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <array>
#include <vector>

int N_sizes[] = {512, 2048, 1024, 1024, 2048, 2048, 512, 1024, 1024, 2048, 1024, 1024, 2048, 512, 512, 1024, 512, 2048, 512, 512, 2048};
int M_sizes[] = {2048, 1024, 512, 512, 512, 1024, 512, 2048, 2048, 1024, 1024, 512, 1024, 1024, 512, 2048, 1024, 2048, 2048, 512, 2048};
int K_sizes[] = {1024, 1024, 512, 1024, 512, 2048, 512, 2048, 1024, 512, 2048, 512, 2048, 1024, 1024, 512, 1024, 512, 2048, 2048, 2048};


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

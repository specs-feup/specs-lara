cmake_minimum_required(VERSION 2.8)
project (mm)

add_executable(mm matrix_mul.cpp matrix_mul.hpp)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -fopenmp")

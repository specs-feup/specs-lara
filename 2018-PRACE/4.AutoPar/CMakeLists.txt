################ standard cmake script ################
cmake_minimum_required(VERSION 3.0)
project(AutoPar C)

set(CMAKE_C_STANDARD 11)

# OpenMP flags
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    
add_executable(matrix_mul src/matrix_mul.c)

################ clava-specific instructions ################
find_package(Clava REQUIRED)

clava_weave(matrix_mul lara/AutoPar.lara)


################ standard cmake script ################
cmake_minimum_required(VERSION 3.0)
project(Exploration C)

set(CMAKE_C_STANDARD 11)

# OpenMP flags
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    
add_executable(matrix_mul src/matrix_mul.c)

################ clava-specific instructions ################
# Put your own code here!

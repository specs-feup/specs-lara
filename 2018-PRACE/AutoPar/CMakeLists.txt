# standard cmake script
cmake_minimum_required(VERSION 3.0)
project(AutoPar C)

set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
add_executable(matrix_mul src/matrix_mul.c)

# clava-specific instructions
find_package(Clava REQUIRED)

clava_weave(matrix_mul AutoPar.lara)

	find_package(OpenMP REQUIRED)
	
	
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

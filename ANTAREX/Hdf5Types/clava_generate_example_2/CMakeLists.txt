cmake_minimum_required(VERSION 3.3)
project(test_lib)

set(Clava_DIR "C:/Users/JoaoBispo/Desktop/shared/repositories-programming/clava/CMake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fopenmp -Wno-invalid-offsetof")

set(SOURCE_FILES
		src/test_lib.cpp
)


add_library(test_lib ${SOURCE_FILES})

	
# Add include folders
target_include_directories(test_lib PUBLIC "include")

find_package(Clava REQUIRED)
#clava_generate(test_lib test_lib_generated "Generate.lara" ARGS "a: 10, b: 'hello'")  
clava_generate_hdf5(test_lib test_lib_hdf5 FILTER "kind: 'class', rx_name: 'A'" "include/test_lib.h")  

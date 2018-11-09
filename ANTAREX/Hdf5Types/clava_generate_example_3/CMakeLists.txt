cmake_minimum_required(VERSION 3.3)
project(test_lib CXX)

set(Clava_DIR "C:/Users/JoaoBispo/Desktop/shared/repositories-programming/clava/CMake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fopenmp -Wno-invalid-offsetof")

add_library(test_headers include/test_lib.h)

# Header-based target, needs language to be set otherwise CMake will not be able to determine it
set_target_properties(test_headers PROPERTIES LINKER_LANGUAGE CXX)

# Add include folders
target_include_directories(test_headers PUBLIC "include")

find_package(Clava REQUIRED)
clava_generate_hdf5(test_headers test_hdf5 FILTER "kind: 'class', rx_name: 'A'")
clava_generate(test_headers test_generated "Generate.lara")  

add_executable(test src/test_lib.cpp)

target_include_directories(test PUBLIC  ${test_hdf5_INCLUDE_DIRS} "include")
target_link_libraries(test ${test_hdf5_LIBRARIES})


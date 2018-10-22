cmake_minimum_required(VERSION 2.8.12)
project(HDF5 CXX)

# Required if CMake plugin is not installed
set(Clava_DIR "C:/Users/JoaoBispo/Desktop/shared/repositories-programming/clava/CMake")

add_library(h5ex_headers Data.h)

# Header-based target, needs language to be set otherwise CMake will not be able to determine it
set_target_properties(h5ex_headers PROPERTIES LINKER_LANGUAGE CXX)

find_package(Clava REQUIRED)
clava_generate_hdf5(h5ex_headers h5ex_hdf5)

add_executable(h5ex main.cpp)
target_link_libraries(h5ex ${h5ex_hdf5_LIBRARIES})
target_include_directories(h5ex PUBLIC  ${h5ex_hdf5_INCLUDE_DIRS})

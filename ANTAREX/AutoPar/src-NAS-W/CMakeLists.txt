cmake_minimum_required(VERSION 3.5)
project(nas_w)

set(Clava_DIR C:/Users/JoaoBispo/Desktop/shared/repositories-programming/clava/CMake)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O3 -mcmodel=medium")

add_subdirectory(common)
#add_subdirectory(BT)
#add_subdirectory(CG)
#add_subdirectory(DC)
add_subdirectory(EP)



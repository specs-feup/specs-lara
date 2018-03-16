cmake_minimum_required(VERSION 3.0.2)

set (EXE_NAME "variant_0")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

set (CMAKER_SOURCES "/home/pedro/Documents/repositories/specs-lara/ANTAREX/AutotunerMatrixMult/date/src/lat_output/matrix_mul/variants/0/src/matrix_mul.hpp"
                    "/home/pedro/Documents/repositories/specs-lara/ANTAREX/AutotunerMatrixMult/date/src/lat_output/matrix_mul/variants/0/src/matrix_mul.cpp"
)

include_directories("../0")

add_executable(${EXE_NAME} ${CMAKER_SOURCES})

target_link_libraries(${EXE_NAME} "stdc++")


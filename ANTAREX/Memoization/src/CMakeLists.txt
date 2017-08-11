cmake_minimum_required(VERSION 3.2)    
    
project(MEMOI_TEST C)

set(CMAKE_C_STANDARD 11)

set(MEMOI_TEST_SRC log_op.c)

add_executable(memoi_test ${MEMOI_TEST_SRC})

target_link_libraries(memoi_test m)

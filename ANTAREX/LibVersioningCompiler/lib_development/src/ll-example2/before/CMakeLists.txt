cmake_minimum_required(VERSION 3.0.2)

set (EXE_NAME "Example2")

set (SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set (SRC_CPP ${SRC_DIR}/main.cpp
             ${SRC_DIR}/kernel.cpp
    )

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set (DEP_LIBS "stdc++")

set (DEP_INC_DIR "")

include_directories(${DEP_INC_DIR})
add_executable(${EXE_NAME} ${SRC_CPP})
target_link_libraries(${EXE_NAME} ${DEP_LIBS})

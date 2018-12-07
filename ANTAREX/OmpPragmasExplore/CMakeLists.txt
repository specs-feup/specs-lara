cmake_minimum_required(VERSION 3.5)
project(gemm)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
#set(Clava_DIR <CMAKE_CLAVA_PLUGIN>)
set(Clava_DIR C:/Users/JoaoBispo/Desktop/shared/repositories-programming/clava/CMake)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O3 -mcmodel=medium")

set(SOURCE_FILES
	src/gemm.c
)


add_executable(gemm ${SOURCE_FILES})

#target_link_libraries(gemm utilities)
#add_subdirectory(utilities)

# OpenMP Support
#find_package(OpenMP REQUIRED)
	
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" )
	
# Clava section

find_package(Clava REQUIRED)

#clava_weave(gemm ReenablePragmas.lara)
clava_weave(gemm ExplorePragmas.lara)
#clava_weave(gemm SelectPragmasV2.lara)


# Test

#add_custom_command(
#     TARGET gemm
#     COMMENT "Running compiled application"
#     POST_BUILD 
#     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#     COMMAND gemm
#)
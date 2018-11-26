cmake_minimum_required(VERSION 3.5)
project(rebuild_lib)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
#set(Clava_DIR "../cmake-clava")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -O3")

set(SOURCE_FILES
		rebuild_lib.c
)



add_library(rebuild_lib ${SOURCE_FILES})

# Add include folders
#target_include_directories(rebuild_lib PUBLIC ../library "./")
target_include_directories(rebuild_lib PUBLIC "./")

add_subdirectory("../library_pre" "${CMAKE_CURRENT_BINARY_DIR}/rebuild_library_pre_build")
target_link_libraries(rebuild_lib rebuild_lib_pre)

# Clava section

#find_package(Clava REQUIRED)

# Use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
#clava_weave_autopar(rebuild_lib)

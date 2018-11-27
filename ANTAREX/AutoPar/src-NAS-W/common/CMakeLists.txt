cmake_minimum_required(VERSION 3.5)
project(nas_common)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
# set(Clava_DIR <PATH TO CLAVA CMAKE PLUGIN FOLDER>)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

set(SOURCE_FILES
	c_print_results.c
	c_timers.c
	print_results.c
	randdp.c
	wtime.c
	
)



add_library(nas_common ${SOURCE_FILES})

	
# Add include folders
target_include_directories(nas_common PUBLIC "./")

# Clava section

#find_package(Clava REQUIRED)

# Use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
#clava_weave_autopar(bt_w)

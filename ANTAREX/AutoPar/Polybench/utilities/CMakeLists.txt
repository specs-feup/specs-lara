cmake_minimum_required(VERSION 3.5)
project(utilities)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
#set(Clava_DIR <CMAKE_CLAVA_PLUGIN>)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O3 -mcmodel=medium")

set(SOURCE_FILES
	polybench.c
)


add_library(utilities ${SOURCE_FILES})

# Add include folders
target_include_directories(utilities PUBLIC "./")	

# Just in case a benchmark needs it
target_link_libraries(utilities m)

cmake_minimum_required(VERSION 3.5)
project(dc_w)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
#set(Clava_DIR <CMAKE_CLAVA_PLUGIN>)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -O3 -mcmodel=medium")

set(SOURCE_FILES
		adc.c
		dc.c
		extbuild.c
		jobcntl.c
		rbt.c
)


add_executable(dc_w ${SOURCE_FILES})

# Add nas_common
#add_subdirectory("../common" "${CMAKE_CURRENT_BINARY_DIR}/nas_common_build")
target_link_libraries(dc_w nas_common m)

# Add include folders
target_include_directories(dc_w PUBLIC "./")

	
# Clava section

find_package(Clava REQUIRED)

# Use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
clava_weave_autopar(dc_w)


# Test

add_custom_command(
     TARGET dc_w
     COMMENT "Running compiled application"
     POST_BUILD 
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
     COMMAND dc_w
)

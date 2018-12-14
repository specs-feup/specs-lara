cmake_minimum_required(VERSION 3.5)
project(ft_w)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
#set(Clava_DIR <CMAKE_CLAVA_PLUGIN>)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O3 -mcmodel=medium")

set(SOURCE_FILES
		appft.c
		auxfnct.c
		fft3d.c
		mainft.c
		verify.c
)


add_executable(ft_w ${SOURCE_FILES})

# Add nas_common
target_link_libraries(ft_w nas_common m)

# Add include folders
target_include_directories(ft_w PUBLIC "./")

	
# Clava section

find_package(Clava REQUIRED)

# Use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
clava_weave_autopar(ft_w)


# Test

add_custom_command(
     TARGET ft_w
     COMMENT "Running compiled application"
     POST_BUILD 
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
     COMMAND ft_w
)
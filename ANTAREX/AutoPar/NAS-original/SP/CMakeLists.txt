cmake_minimum_required(VERSION 3.5)
project(sp_w)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
#set(Clava_DIR <CMAKE_CLAVA_PLUGIN>)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O3 -mcmodel=medium")

set(SOURCE_FILES
		add.c
		adi.c
		error.c
		exact_rhs.c
		exact_solution.c
		initialize.c	
		ninvr.c		
		pinvr.c
		rhs.c
		set_constants.c
		sp.c
		txinvr.c
		tzetar.c
		verify.c
		x_solve.c
		y_solve.c
		z_solve.c
)


add_executable(sp_w ${SOURCE_FILES})

# Add nas_common
target_link_libraries(sp_w nas_common m)


# Add include folders
target_include_directories(sp_w PUBLIC "./")

	
# Clava section

find_package(Clava REQUIRED)

# Use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
clava_weave_autopar(sp_w)


# Test

add_custom_command(
     TARGET sp_w
     COMMENT "Running compiled application"
     POST_BUILD 
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
     COMMAND sp_w
)
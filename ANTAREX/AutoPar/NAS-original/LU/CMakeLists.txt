cmake_minimum_required(VERSION 3.5)
project(lu_w)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
#set(Clava_DIR <CMAKE_CLAVA_PLUGIN>)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -O3 -mcmodel=medium")

set(SOURCE_FILES
		blts.c
		buts.c
		domain.c
		erhs.c
		error.c
		exact.c
		jacld.c
		jacu.c
		l2norm.c
		lu.c
		pintgr.c
		read_input.c
		rhs.c
		setbv.c
		setcoeff.c
		setiv.c
		ssor.c
		verify.c
)


add_executable(lu_w ${SOURCE_FILES})

# Add nas_common
target_link_libraries(lu_w nas_common m)


# Add include folders
target_include_directories(lu_w PUBLIC "./")

	
# Clava section

find_package(Clava REQUIRED)

# Use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
clava_weave_autopar(lu_w)


# Test

add_custom_command(
     TARGET lu_w
     COMMENT "Running compiled application"
     POST_BUILD 
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
     COMMAND lu_w
)
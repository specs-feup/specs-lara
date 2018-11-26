cmake_minimum_required(VERSION 3.5)
project(bt_w)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
#set(Clava_DIR <CMAKE_CLAVA_PLUGIN>)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -O3 -mcmodel=medium")

set(SOURCE_FILES
		add.c
		adi.c
		bt.c
		error.c
		exact_rhs.c
		exact_solution.c
		initialize.c		
		rhs.c
		set_constants.c
		solve_subs.c
		verify.c
		x_solve.c
		y_solve.c
		z_solve.c
		
		# Add headers, for auto-parallelization to have access to the code
		#header.h
		#npbparams.h
		#work_lhs.h
)



add_executable(bt_w ${SOURCE_FILES})

add_subdirectory("../common" "${CMAKE_CURRENT_BINARY_DIR}/nas_common_build")
target_link_libraries(bt_w nas_common m)
#target_include_directories(bt_w PUBLIC ../common "./")
target_include_directories(bt_w PUBLIC "./")

	
# Add include folders
#target_include_directories(bt_w PUBLIC "./")

# Clava section

find_package(Clava REQUIRED)

# Use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
clava_weave_autopar(bt_w)

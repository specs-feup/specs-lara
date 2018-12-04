cmake_minimum_required(VERSION 3.5)
project(autopar_test)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
# set(Clava_DIR <PATH TO CLAVA CMAKE PLUGIN FOLDER>)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-invalid-offsetof")

set(SOURCE_FILES
		src/main.cpp
)

add_executable(autopar_test ${SOURCE_FILES})

	
# Add include folders
target_include_directories(autopar_test PUBLIC "include")

# Clava section

find_package(Clava REQUIRED)

# Apply a generic LARA file to the current code, transforming it
clava_weave(autopar_test "ApplyAutopar.lara")

# Alternatively, use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
#clava_weave_autopar(autopar_test)

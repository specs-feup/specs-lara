cmake_minimum_required(VERSION 3.3)
project(autopar_test)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
#set(Clava_DIR "<PATH_TO_CLAVA_CMAKE_PLUGIN_FOLDER>")

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
#clava_weave(autopar_test "ApplyAutopar.lara")

# Alternatively, use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
clava_weave_autopar(autopar_test)

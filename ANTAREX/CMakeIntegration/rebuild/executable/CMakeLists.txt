cmake_minimum_required(VERSION 3.5)
project(rebuild_exe)

# If Clava CMake plugin is not installed (e.g., sudo clava-update), the path to the plugin
# can be manually specified by setting the variable Clava_DIR
# The Clava CMake plugin can be found here: https://github.com/specs-feup/clava/tree/master/CMake
#set(Clava_DIR <CLAVA CMAKE DIR>)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -O3")

set(SOURCE_FILES
		rebuild.c
)

add_executable(rebuild_exe ${SOURCE_FILES})

add_subdirectory("../library" "${CMAKE_CURRENT_BINARY_DIR}/rebuild_library_build")
target_link_libraries(rebuild_exe rebuild_lib m)
 
# Add include folders
target_include_directories(rebuild_exe PUBLIC "./")

# Clava section

find_package(Clava REQUIRED)

clava_weave(rebuild_exe Rebuild.lara)


# Test

add_custom_command(
     TARGET rebuild_exe
     COMMENT "Running compiled application"
     POST_BUILD 
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
     COMMAND rebuild_exe
)
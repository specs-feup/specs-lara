# check the minimum version
cmake_minimum_required( VERSION 3.0.2 )

# the project name
project( mmul )

################################
#### General configure section
################################

# force the Release build if not already set
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif(NOT CMAKE_BUILD_TYPE)

# setting debug flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g3 -O0")

# setting release with debug info flags
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -mtune=native -g3 -O2")

# setting release flags
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native -O3")



################################
#### External dipendency handler
################################

# add the path to the modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/high_level_interface/cmake")

# find the required packages
find_package(MARGOT_HEEL)

# add the required include directories
include_directories(${MARGOT_HEEL_INCLUDES})



################################
#### Build section
################################


# Set the application name
set( APP_NAME "mm")

# build the application
add_executable(${APP_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/matrix_mult/matmult.c)

# make sure to link with g++
set_target_properties(${APP_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# link the high level interface
target_link_libraries(${APP_NAME} ${MARGOT_HEEL_LIBRARIES})

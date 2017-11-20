# check the minimum version
cmake_minimum_required( VERSION 2.8.12 )

# the project name
project( MATRIX_MUL )

################################
#### General configure section
################################

# force the Release build if not already set
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# setting common c++ flags
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -march=native -fopenmp" )

# setting debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g3 -O3 -pg")

# setting release with debug info flags
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -mtune=native -g3 -O2")

# setting release flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native -O3")

# disable some linux warning
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result -Wno-format" )



# Set the application name
set( APP_NAME "matrix_mul")

#set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/weaved_code)
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

#set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/weaved_code/include)
#set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add the include folder
#include_directories( ${INCLUDE_DIR})


# Source files
set(SRC_FILES
	${SOURCE_DIR}/matrix_mul.c
)

# Libraries
set(EXTERNAL_LIB_LINK_LIB )

add_executable(${APP_NAME} ${SRC_FILES})
target_link_libraries(${APP_NAME} ${EXTERNAL_LIB_LINK_LIB} )

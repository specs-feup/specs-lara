cmake_minimum_required(VERSION 3.0)

# Every project needs a name.
project(aos2soa)

option(LVL "lvl" 0)
option(AFTER "after" OFF)


message(STATUS "[AoS2SoA] Setting source level to ${LVL}")
message(STATUS "[AoS2SoA] Use 'after' version: ${AFTER}")


set(CMAKE_CXX_FLAGS "-O3")

set(SOURCE_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/src-lvl${LVL}" )
if(AFTER)
	set(SOURCE_FOLDER "${SOURCE_FOLDER}-after" )
endif()

include_directories(${SOURCE_FOLDER})
include_directories("include")


# set the sources of the application
set( SOURCE
	${SOURCE_FOLDER}/aos2soa.cpp
	${SOURCE_FOLDER}/main.cpp
)

# compile the application
add_executable("aos2soa" ${SOURCE})

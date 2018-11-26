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
#target_link_libraries(rebuild_exe rebuild_lib)
 
# Add include folders
#target_include_directories(rebuild_exe PUBLIC ../library "./")
target_include_directories(rebuild_exe PUBLIC "./")

#get_target_property(INC rebuild_exe INCLUDE_DIRECTORIES)
#get_target_property(INTERFACE_INC rebuild_exe INTERFACE_INCLUDE_DIRECTORIES)
#get_target_property(INTERFACE_SYSTE_INC rebuild_exe INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
#get_target_property(C_INC rebuild_exe C_INCLUDE_WHAT_YOU_USE)
#get_target_property(CXX_INC rebuild_exe CXX_INCLUDE_WHAT_YOU_USE)

#message(STATUS "INC: ${INC}")
#message(STATUS "INTERFACE_INC: ${INTERFACE_INC}")
#message(STATUS "INTERFACE_SYSTE_INC: ${INTERFACE_SYSTE_INC}")
#message(STATUS "C_INC: ${C_INC}")
#message(STATUS "CXX_INC: ${CXX_INC}")


#get_target_property(LINK rebuild_exe LINK_LIBRARIES)
#get_target_property(INTERFACE_LINK rebuild_exe LINK_INTERFACE_LIBRARIES)

#message(STATUS "LINK: ${LINK}")
#message(STATUS "INTERFACE_LINK: ${INTERFACE_LINK}")
#if("${LINK}" MATCHES "LINK-NOTFOUND")
#	set(LINK "")
#endif()

#list(LENGTH LINK LINK_LENGTH)
#message(STATUS "LINK LENGTH: ${LINK_LENGTH}")
	

#set(ORIG_INCLUDES "./")	
	
#foreach(ORIG_LINK_LIB IN LISTS LINK)
#		#message(STATUS "ORIG_LINK_LIB: ${ORIG_LINK_LIB}")
#		get_target_property(ORIG_LINK_LIB_INC ${ORIG_LINK_LIB} INCLUDE_DIRECTORIES)
#		if((${ORIG_LINK_LIB_INC} MATCHES "ORIG_LINK_LIB_INC-NOTFOUND"))
#			set(ORIG_LINK_LIB_INC "")
#		endif()
#		list(APPEND ORIG_INCLUDES ${ORIG_LINK_LIB_INC})
#		message(STATUS "ORIG_LINK_LIB_INC: ${ORIG_LINK_LIB_INC}")
#endforeach(ORIG_LINK_LIB)
#message(STATUS "ORIG_INCLUDES: ${ORIG_INCLUDES}")

		
# Clava section

find_package(Clava REQUIRED)

#set(INCLUDES_RECURSIVE "")
#target_include_directories_recursive(rebuild_exe INCLUDES_RECURSIVE)
#message(STATUS "INCLUDES_RECURSIVE: ${INCLUDES_RECURSIVE}")

clava_weave(rebuild_exe Rebuild.lara)

# Use the utility call for AutoPar.
# This function automatically adds the flags required for compiling OpenMP
#clava_weave_autopar(rebuild_exe)

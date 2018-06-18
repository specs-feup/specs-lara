cmake_minimum_required (VERSION 3.0)

project (synt)

find_package(Clava REQUIRED)

# TO MOVE TO LIBRARY
function(clava_weave ORIG_TARGET ASPECT)

	# get CMakeLists.txt dir
	get_target_property(ORIG_CMAKE_DIR ${ORIG_TARGET} SOURCE_DIR)
	message(STATUS "ORIG_CMAKE_DIR: ${ORIG_CMAKE_DIR}")
	
	# get full path of aspect file
	set(ASPECT "${ORIG_CMAKE_DIR}/${ASPECT}")
	message(STATUS "ASPECT: ${ASPECT}")

	# get woven directory path
	set(WOVEN_DIR "clava_${ORIG_TARGET}")
	message(STATUS "WOVEN_DIR: ${WOVEN_DIR}")

	# get original source files
	get_target_property(ORIG_SOURCES ${ORIG_TARGET} SOURCES)
	message(STATUS "ORIG_SOURCES: ${ORIG_SOURCES}")
	
	# process original source file list
	string(REGEX REPLACE "([^;]+)" "${ORIG_CMAKE_DIR}/\\1" PROC_ORIG_SOURCES "${ORIG_SOURCES}")
	if(NOT "${CMAKE_HOST_SYSTEM}" MATCHES ".*Windows.*")
		string(REGEX REPLACE ";" ":" PROC_ORIG_SOURCES "${PROC_ORIG_SOURCES}")
	endif()
	message(STATUS "PROC_ORIG_SOURCES: ${PROC_ORIG_SOURCES}")
	
	# get original include directories
	get_target_property(ORIG_INC_DIRS ${ORIG_TARGET} INCLUDE_DIRECTORIES)
	message(STATUS "ORIG_INC_DIRS: ${ORIG_INC_DIRS}")
	
	# make the cmake configuration depend on the LARA file
	set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${ASPECT};${ORIG_SOURCES}")
	
	# execute Clava (TODO: set correct standard from cmake) (TODO: set correct clava jar)
	execute_process(
		COMMAND clava ${ASPECT} --cmake -b 2 -std c99 -p "${PROC_ORIG_SOURCES}" -of "${WOVEN_DIR}"
		WORKING_DIRECTORY ${ORIG_CMAKE_DIR}
	)
	
	# read new sources
	if(EXISTS "${ORIG_CMAKE_DIR}/${WOVEN_DIR}/clava_generated_files.txt")
        file(READ "${ORIG_CMAKE_DIR}/${WOVEN_DIR}/clava_generated_files.txt" CLAVA_WOVEN_SOURCES)
        string(STRIP "${CLAVA_WOVEN_SOURCES}" CLAVA_WOVEN_SOURCES)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CLAVA_WOVEN_SOURCES}")

        message(STATUS "CLAVA_WOVEN_SOURCES: ${CLAVA_WOVEN_SOURCES}")
	else()
		# TODO: add something here
		message(STATUS "couldn't find sources")
	endif()

	# set new sources
	set_target_properties(${ORIG_TARGET}
		PROPERTIES SOURCES "${CLAVA_WOVEN_SOURCES}"
	)
	
endfunction(clava_weave)
# TO MOVE TO LIBRARY


add_executable(synt src/main.c src/print.c src/print.h)

clava_weave(synt "lara/Synt.lara")

cmake_minimum_required(VERSION 3.3)
project(HDF5_Clava)

find_package(HDF5 COMPONENTS "CXX" REQUIRED)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
endif()

# Enums headers
set(FOLDER_ENUMS src/Enums)
include_directories(${FOLDER_ENUMS})
set(SOURCE_FILES_ENUMS
        ${FOLDER_ENUMS}/Direction.h
		${FOLDER_ENUMS}/ESpecificInfo.h
		${FOLDER_ENUMS}/Direction.h)

# Structures headers
set(FOLDER_STRUCTURES src/Structures)
include_directories(${FOLDER_STRUCTURES})
set(SOURCE_FILES_STRUCTURES
    ${FOLDER_STRUCTURES}/EdgeData.h
		${FOLDER_STRUCTURES}/EdgeRaw.h
		${FOLDER_STRUCTURES}/GeometryShortcut.h
		${FOLDER_STRUCTURES}/IndexGeometryShortcut.h
		${FOLDER_STRUCTURES}/NodePosition.h
		${FOLDER_STRUCTURES}/NodeRaw.h
		${FOLDER_STRUCTURES}/RestrictedTurn.h
		${FOLDER_STRUCTURES}/RoutingIndexInfo.h
		${FOLDER_STRUCTURES}/ShortcutPosition.h
		${FOLDER_STRUCTURES}/Geo/GpsPoint.h)

# Utility headers
set(FOLDER_UTILITY src/Utility)
include_directories(${FOLDER_UTILITY})
set(SOURCE_FILES_UTILITY
        ${FOLDER_UTILITY}/Definitions.h)

set(HEADER_FILES
		${SOURCE_FILES_ENUMS}
		${SOURCE_FILES_STRUCTURES}
		${SOURCE_FILES_UTILITY})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fopenmp -Wno-invalid-offsetof")

set(SOURCE_FILES
		src-test/test.cpp
		${HEADER_FILES})

#message(STATUS "BIN -> ${CMAKE_CURRENT_BINARY_DIR}")
#message(STATUS "SOURCE -> ${CMAKE_SOURCE_DIR}")
#message(STATUS "CURRENT SOURCE -> ${CMAKE_CURRENT_SOURCE_DIR}")

set(CLAVA_COMMAND
		#clava -c ${CMAKE_SOURCE_DIR}/hdf5.clava -av "\'${CMAKE_SOURCE_DIR}\'"
		clava -c ${CMAKE_SOURCE_DIR}/hdf5.clava -av \"{srcFolder:'${CMAKE_SOURCE_DIR}'}\"
)

message(STATUS "CLAVA Command: ${CLAVA_COMMAND}")

add_custom_command(
	OUTPUT ${CMAKE_SOURCE_DIR}/generated/CompType.cpp ${CMAKE_SOURCE_DIR}/generated/CompType.h
	#COMMAND clava -c ${CMAKE_SOURCE_DIR}/hdf5.clava -av "\'${CMAKE_SOURCE_DIR}\'"
	COMMAND ${CLAVA_COMMAND}
	COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/weaved_code/lara-generated/CompType.cpp ${CMAKE_SOURCE_DIR}/generated/CompType.cpp
	COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/weaved_code/lara-generated/CompType.h ${CMAKE_SOURCE_DIR}/generated/CompType.h
	COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_SOURCE_DIR}/generated
	DEPENDS ${HEADER_FILES}
)

add_custom_target(hdf5_converters ALL DEPENDS ${CMAKE_SOURCE_DIR}/generated/CompType.cpp ${CMAKE_SOURCE_DIR}/generated/CompType.h)

add_executable(hdf5_test
	${SOURCE_FILES}
	${CMAKE_SOURCE_DIR}/generated/CompType.cpp
	${CMAKE_SOURCE_DIR}/generated/CompType.h)

add_dependencies(hdf5_test hdf5_converters)
target_include_directories(hdf5_test PUBLIC ${CMAKE_SOURCE_DIR}/generated)
target_link_libraries(hdf5_test ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})

cmake_minimum_required(VERSION 3.3)
project(test_lib CXX)

set(Clava_DIR "C:/Users/JoaoBispo/Desktop/shared/repositories-programming/clava/CMake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fopenmp -Wno-invalid-offsetof")


add_library(test_headers include/test_lib.h)

# Header-based target, needs language to be set otherwise CMake will not be able to determine it
set_target_properties(test_headers PROPERTIES LINKER_LANGUAGE CXX)

#get_target_property(TEST_HEADERS_INCLUDES_BEFORE test_headers INCLUDE_DIRECTORIES)

# Add include folders
target_include_directories(test_headers PUBLIC "include")

#get_target_property(TEST_HEADERS_INCLUDES_AFTER test_headers INCLUDE_DIRECTORIES)
#message(STATUS "TEST_HEADER Includes before: ${TEST_HEADERS_INCLUDES_BEFORE}")
#message(STATUS "TEST_HEADER Includes after: ${TEST_HEADERS_INCLUDES_AFTER}")

find_package(Clava REQUIRED)
#clava_generate(test_lib test_lib_generated "Generate.lara" "a: 10, b: 'hello'")  
clava_generate_hdf5(test_headers test_hdf5 "kind: 'class', rx_name: 'A'")  

clava_generate(test_headers test_generated "Generate.lara")  

# Add HDF5 dependency to test_hdf5
#find_package(HDF5 REQUIRED)

# Generate HDF5 interfaces for current code
#target_include_directories(test_hdf5 PUBLIC ${HDF5_INCLUDE_DIRS})
#target_link_libraries(test_hdf5 ${HDF5_LIBRARIES})
# Should use HDF5 variable with link flags
#target_link_libraries(test_hdf5 hdf5_cpp)



add_executable(test src/test_lib.cpp)

#get_target_property(TEST_INCLUDES_BEFORE test INCLUDE_DIRECTORIES)

#add_dependencies(test test_hdf5)
#get_target_property(TEST_HEADERS_INCLUDES_AFTER test_headers INCLUDE_DIRECTORIES)
#get_target_property(TEST_HEADERS_INCLUDES test_headers INCLUDE_DIRECTORIES)


#target_include_directories(test PUBLIC ${TEST_HEADERS_INCLUDES} ${TEST_HDF5_INCLUDES})	




get_target_property(TEST_HDF5_INCLUDES test_hdf5 INCLUDE_DIRECTORIES)
target_include_directories(test PUBLIC ${TEST_HDF5_INCLUDES})

#target_include_directories(test PUBLIC test_hdf5)

#get_target_property(TEST_INCLUDES test INCLUDE_DIRECTORIES)
#message(STATUS "TEST INCLUDES:  ${TEST_INCLUDES}")

get_target_property(TEST_HDF5_LINK_LIBRARIES test_hdf5 LINK_LIBRARIES)
target_link_libraries(test ${TEST_HDF5_LINK_LIBRARIES})
#target_link_libraries(test test_hdf5)

#add_dependencies(test test_hdf5)




#target_link_libraries(test ${TEST_HDF5_LINK_LIBRARIES})
#target_link_libraries(test ${HDF5_LIBRARIES})
#target_link_libraries(test hdf5_cpp)
#message(STATUS "HDF5 LIBS: ${HDF5_LIBRARIES}")
#get_target_property(TEST_HEADERS_INCLUDES_AFTER_2 test_headers INCLUDE_DIRECTORIES)
#message(STATUS "TEST_HEADER Includes after 2: ${TEST_HEADERS_INCLUDES_AFTER_2}")


#get_target_property(TEST_INCLUDES_AFTER test INCLUDE_DIRECTORIES)
#message(STATUS "TEST Includes before: ${TEST_INCLUDES_BEFORE}")
#message(STATUS "TEST Includes after: ${TEST_INCLUDES_AFTER}")
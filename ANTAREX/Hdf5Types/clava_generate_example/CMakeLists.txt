cmake_minimum_required(VERSION 3.3)
project(test_lib)

set(Clava_DIR "C:/Users/JoaoBispo/Desktop/shared/repositories-programming/clava/CMake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fopenmp -Wno-invalid-offsetof")

set(SOURCE_FILES
		src/test_lib.cpp
)


add_library(test_lib ${SOURCE_FILES})

	
# Add include folders
target_include_directories(test_lib PUBLIC "include")

find_package(Clava REQUIRED)
clava_generate(test_lib test_lib_generated ARGS "Generate.lara" "a: 10, b: 'hello'")  
#clava_generate_hdf5(test_lib test_lib_hdf5 FILTER "rx_name: 'A', kind: 'class'")

#clava_generate_hdf5(test_lib test_lib_hdf5)
#clava_generate_hdf5(test_lib test_lib_hdf5 FILTER "rx_name: 'A'")

# HDF5 Support
#find_package(HDF5 REQUIRED)
#find_package(Clava REQUIRED)

# Generate HDF5 interfaces for current code
#clava_generate(test_lib test_lib_hdf5 "ApplyHdf5.lara" ARGS "filter: {rx_name: 'A', kind: 'class'}") # If attribute is prefixed by 'rx_', it is interpreted as a regex 
#target_include_directories(test_lib_hdf5 PUBLIC ${HDF5_INCLUDE_DIRS})
#target_link_libraries(test_lib_hdf5 ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})


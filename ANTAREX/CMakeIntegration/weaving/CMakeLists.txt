cmake_minimum_required(VERSION 3.3)
project(test)

set(Clava_DIR "C:/Users/JoaoBispo/Desktop/shared/repositories-programming/clava/CMake")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")

set(SOURCE_FILES main.c)

add_executable(test ${SOURCE_FILES})

find_package(Clava REQUIRED)

# Apply a generic LARA file to the current code, transforming it
#clava_weave(test "Test.lara" "input1: '1st insert'")
#clava_weave(test "Test.lara" "input1: '2nd insert'")
clava_weave(test "Test.lara" ARGS "input1: '1st insert'" FLAGS -b 3)
clava_weave(test "Test.lara" ARGS "input1: '2nd insert'")


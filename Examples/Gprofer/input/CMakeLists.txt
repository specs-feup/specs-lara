cmake_minimum_required(VERSION 2.8)

set (EXE_NAME "gprofer_bin")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-pie -pg")

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -no-pie -pg")

set (CMAKER_SOURCES "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/disparity_matisse_modified.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/getDisparity.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/lib/array_creators_alloc.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/lib/general_matrix.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/lib/load.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/lib/matisse_string.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/lib/matlab_general.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/lib/matrix.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/lib/tensor.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/main_test.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/new_from_matrix_double_0_d.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/padarray_specialized_pre_post_tdd1x2_2.c"
                    "/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src/padarray_specialized_pre_tdi1x2s_2.c"
)

include_directories("/home/pedro/Documents/repositories/specs-lara/Examples/Gprofer/input/src")

add_executable(${EXE_NAME} ${CMAKER_SOURCES})

target_link_libraries(${EXE_NAME} "")


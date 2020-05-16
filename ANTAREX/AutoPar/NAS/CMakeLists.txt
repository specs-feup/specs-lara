cmake_minimum_required(VERSION 3.5)
project(nas_single_source_w)

# Putting a smaller example first
add_subdirectory(EP)

add_subdirectory(BT)
add_subdirectory(CG)
add_subdirectory(FT)
add_subdirectory(IS)
add_subdirectory(LU)
add_subdirectory(MG)
add_subdirectory(SP)
add_subdirectory(UA)

cmake_minimum_required(VERSION 3.5)
project(ANTAREX-TESTS)

# AutoPar
#add_subdirectory(AutoPar/NAS)
add_subdirectory(AutoPar/Polybench)

# CMake Plugin
add_subdirectory(CMakeIntegration/rebuild/executable)


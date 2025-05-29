cmake_minimum_required(VERSION 3.10)
project(ANTAREX-TESTS)

# AutoPar tests
add_subdirectory(AutoPar)

# CMake Plugin tests
add_subdirectory(CMakeIntegration/rebuild/executable)


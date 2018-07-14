cmake_minimum_required (VERSION 3.0)

project (synt)

find_package(Clava REQUIRED)

set(SOURCES src/main.c src/print.c src/print.h)

add_executable(synt "${SOURCES}")

clava_weave(synt "lara/Synt.lara")

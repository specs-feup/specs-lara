cmake_minimum_required (VERSION 3.0)

project (synt)

find_package(Clava REQUIRED)

add_executable(synt src/main.c src/print.c src/print.h)

clava_weave(synt "lara/Synt.lara")

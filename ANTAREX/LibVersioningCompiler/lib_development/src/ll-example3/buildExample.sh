#!/bin/bash
cd before
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX="."
make

cd ../..
cd after
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX="."
make

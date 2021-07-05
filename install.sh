#!/bin/bash
if [ ! -d "./build" ]; then
  mkdir ./build
fi
cd ./build

if [ ! -d "build" ]; then
  mkdir build
fi
cd build

# export EIGEN_DIR=../../3rd
cmake ../.. -DEIGEN3_INCLUDE_DIRS=../../3rd
make

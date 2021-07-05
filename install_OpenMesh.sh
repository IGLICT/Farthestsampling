#!/bin/bash
if [ ! -d "./build" ]; then
  mkdir ./build
fi
cd ./build

if [ ! -d "OpenMesh-build" ]; then
  mkdir OpenMesh-build
fi
cd OpenMesh-build

cmake ../../3rd/OpenMesh -DCMAKE_INSTALL_PREFIX=$(pwd)/install-custom -DBUILD_APPS=false
make -j8
make install
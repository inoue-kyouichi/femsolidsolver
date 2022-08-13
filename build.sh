#!/bin/sh
mkdir build
cd build
cmake -D compiler=intel \
      -D CMAKE_INSTALL_PREFIX=/home/totani/bin/femSolidAnalysis \
      -D TP_DIR=/home/totani/lib/TextParser-1.8.5 \
      -D EIGEN_DIR=/home/totani/lib/eigen-3.3.4 \
      -D enable_GLOG=ON \
      -D GLOG_DIR=/home/totani/lib/glog \
      ..

make && make install
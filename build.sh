#!/bin/sh
mkdir build
cd build
cmake -D compiler=intel \
      -D CMAKE_INSTALL_PREFIX=/home/kyouichi/bin/femSolidAnalysis \
      -D TP_DIR=/home/kyouichi/lib/TextParser \
      -D EIGEN_DIR=/home/kyouichi/lib/eigen-3.3.4 \
      -D enable_GLOG=ON \
      -D GLOG_DIR=/home/kyouichi/lib/glog \
      ..

make && make install
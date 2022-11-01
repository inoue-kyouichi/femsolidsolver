#!/bin/sh
mkdir build
cd build
cmake -D compiler=intel \
      -D CMAKE_INSTALL_PREFIX=/home/inoue/bin/femSolidAnalysis \
      -D TP_DIR=/home/inoue/lib/TextParser \
      -D EIGEN_DIR=/home/inoue/lib/eigen-3.3.4 \
      -D enable_GLOG=ON \
      -D GLOG_DIR=/home/inoue/lib/glog \
      ..

make && make install
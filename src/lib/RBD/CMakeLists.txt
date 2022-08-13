###################################################################################
#
# Coherent point drift
#
# Copyright (c) 2020 Biomechanics Lab.,
#                    Department of Mechanical Science and Bioengineering,
#                    Graduate School of Engineering Science,
#                    Osaka University.
# All rights reserved.
#
###################################################################################

set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

set(LIBNAME "RBD")

add_library(${LIBNAME} STATIC)
target_sources(${LIBNAME}
  PRIVATE
    rigidBody.cpp
    TriangleSet.cpp
)

target_include_directories(${LIBNAME} PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/FEM
  ${PROJECT_SOURCE_DIR}/src/lib/linearSolver
)

target_compile_options(${LIBNAME} PUBLIC ${compile_options_tmp})
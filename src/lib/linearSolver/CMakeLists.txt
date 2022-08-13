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

add_library(LS STATIC)
target_sources(LS
  PRIVATE
    pardiso_solver.cpp
    pardiso_solver_rigidBodyInteraction.cpp
)

target_include_directories(LS PRIVATE
  ${PROJECT_SOURCE_DIR}/src/lib/base
)

target_compile_options(LS PUBLIC ${compile_options_tmp})
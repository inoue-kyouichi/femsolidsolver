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

add_library(Base STATIC)
target_sources(Base
  PRIVATE
    math_tools.cpp
    fileIO.cpp
)

target_compile_options(Base PUBLIC ${compile_options_tmp})
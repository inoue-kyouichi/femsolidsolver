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

add_library(FEM STATIC)
target_sources(FEM
  PRIVATE
    domain.cpp
    fem.cpp
    fem_base.cpp
    fem_boundary.cpp
    fem_preprocessing.cpp
    fem_postprocessing.cpp
    fem_SantVenant_spatialForm.cpp
    triangleSurfaceCurvature.cpp
)

target_include_directories(FEM PUBLIC
  ${PROJECT_SOURCE_DIR}/src/lib/base
  ${PROJECT_SOURCE_DIR}/src/lib/linearSolver
)

target_compile_options(FEM PUBLIC ${compile_options_tmp})
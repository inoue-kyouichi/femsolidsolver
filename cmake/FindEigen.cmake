find_path(EIGEN_INCLUDE_DIR Eigen/Core
  PATHS
    ENV EIGEN_ROOT
    ENV EIGEN_INCLUDE_DIR
    ${EIGEN_DIR}
    /usr
    /usr/include/eigen3
    /usr/local
  PATH_SUFFIXES
    include
  )
mark_as_advanced(EIGEN_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen
  REQUIRED_VARS EIGEN_INCLUDE_DIR
  )

if(Eigen_FOUND AND NOT TARGET Eigen::Eigen)
  # EigenはヘッダーのみのライブラリなのでINTERFACEキーワードを指定する
  add_library(Eigen::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN_INCLUDE_DIR}"
    )
  # set_target_propertiesのかわりにtarget_include_directoriesを使ってもOK
  # target_include_directories(Eigen::Eigen INTERFACE ${EIGEN_INCLUDE_DIR})
endif()
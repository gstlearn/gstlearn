project(gstlearn CXX)

cmake_minimum_required(VERSION 3.10)

include(GNUInstallDirs)

set(CMAKE_BUILD_TYPE
  "Release"
  CACHE
  STRING
  "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
  FORCE
)

add_subdirectory(src)
add_subdirectory(doxygen)

install(
  EXPORT corelibs
  FILE gstlearnConfig.cmake
  NAMESPACE gstlearn::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/cmake
)
# Uncomme the following line if you to register the package globally this will
# avoid the manual specification of gstlearn_DIR
# set(CMAKE_EXPORT_PACKAGE_REGISTRY ON) export(PACKAGE gstlearn)

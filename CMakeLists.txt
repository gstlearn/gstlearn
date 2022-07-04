# For recent UseSWIG
cmake_minimum_required(VERSION 3.20)

# Permit using project VERSION
project(gstlearn
        VERSION      0.1.27 # Make it consistent with CHANGES
        DESCRIPTION  "Geostatistics & Machine Learning toolbox"
        HOMEPAGE_URL https://github.com/gstlearn/gstlearn
        LANGUAGES    C CXX) # Enable C language for third party libraries

####################################################
## CONFIGURATION

# By default, activate HDF5 support
option(USE_HDF5 "Using HDF5 support" ON)

# TODO : USE_NETCDF in CMake (mimicking USE_HDF5)
# By default, do not activate NetCDF support
#option(USE_NETCDF "Using NetCDF support" OFF)

# Create gstlearn_DATE, gstlearn_YEAR and gstlearn_MONTH variables
string(TIMESTAMP ${PROJECT_NAME}_DATE "%Y/%m/%d - %H:%M")
string(TIMESTAMP ${PROJECT_NAME}_YEAR "%Y")
string(TIMESTAMP ${PROJECT_NAME}_MONTH "%B")

# Convert project name to uppercase
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UP)

# Fix version.h automatically
configure_file(version.h.in version.h)

# Activate verbose mode
#set(CMAKE_VERBOSE_MAKEFILE ON)

# Options
option(BUILD_PYTHON "Build Python package" ON)
option(BUILD_R      "Build R package"      ON)
if (MINGW)
  set(BUILD_PYTHON OFF)
endif()
if (MSVC)
  set(BUILD_R OFF)
endif()
message(STATUS "BUILD_PYTHON=" ${BUILD_PYTHON})
message(STATUS "BUILD_R="      ${BUILD_R})

####################################################
## BUILD & INSTALL

include(cmake/cpp.cmake)
include(cmake/doxygen.cmake)
include(cmake/doc.cmake)
include(cmake/install.cmake)

####################################################
## PACKAGES

add_subdirectory(python)
add_subdirectory(r)

####################################################
## TESTING

# Add non-regression test target
include(CTest)
enable_testing()

add_subdirectory(tests)
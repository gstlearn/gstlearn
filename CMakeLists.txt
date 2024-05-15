# For recent UseSWIG
cmake_minimum_required(VERSION 3.20)

# Define project here
project(gstlearn
        VERSION      1.2.1 # Make it consistent with CHANGES
        DESCRIPTION  "Geostatistics & Machine Learning toolbox"
        HOMEPAGE_URL "https://gstlearn.org"
        LANGUAGES    C CXX) # Enable C language for third party libraries

####################################################
## CONFIGURATION

# By default, do not activate HDF5 support
option(USE_HDF5 "Using HDF5 support" OFF)

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

# Detect presence of multi configuration generators
get_property(IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

# Activate verbose mode
#set(CMAKE_VERBOSE_MAKEFILE ON)
# Debug find package instruction
#set(CMAKE_FIND_DEBUG_MODE TRUE)

# Options
option(BUILD_PYTHON  "Build Python package"        OFF)
option(BUILD_R       "Build R package"             OFF)
option(BUILD_DOXYGEN "Build Doxygen documentation" OFF)
if (MINGW)
  set(BUILD_PYTHON OFF)
  # TODO : Reactivate HDF5 support for MinGW (currently: link error against szip)
  set(USE_HDF5 OFF)
endif()
if (MSVC)
  set(BUILD_R OFF)
  # TODO : Reactivate HDF5 support for MSVC (currently: hdf5::sz-static target is not found)
  set(USE_HDF5 OFF)
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # Clang compiler detected (define the shortcut variable)
  set(CLANG ON)
endif()
message(STATUS "BUILD_PYTHON="  ${BUILD_PYTHON})
message(STATUS "BUILD_R="       ${BUILD_R})
message(STATUS "BUILD_DOXYGEN=" ${BUILD_DOXYGEN})
message(STATUS "USE_HDF5="      ${USE_HDF5})

# Create gstlearn cmake file path in the build tree
cmake_path(APPEND GSTLEARN_CMAKE_FILE ${PROJECT_BINARY_DIR} "gstlearn.cmake")

####################################################
## EXTERNAL LIBS
add_subdirectory(3rd-party/gmtsph)
add_subdirectory(3rd-party/csparse)

####################################################
## BUILD & INSTALL

include(cmake/version.cmake)
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

####################################################
## DEMONSTRATION SCRIPTS

add_subdirectory(doc/demo)
add_subdirectory(doc/courses)

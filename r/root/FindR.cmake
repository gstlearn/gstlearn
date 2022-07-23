# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# i.e. https://github.com/root-project/root/blob/master/LICENSE
# For the list of contributors see $ROOTSYS/README/CREDITS.

# CMake module to find R
# - Try to find R
# Once done, this will define
#
#  R_FOUND - system has R
#  R_INCLUDE_DIRS - the R include directories
#  R_LIBRARIES - link these to use R
#  R_ROOT_DIR - As reported by R
# Authors: Omar Andres Zapata Mesa 31/05/2013

# Modified by Fabien Ors - 2022
#
# - Retrieve R version and display status message

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_FIND_APPBUNDLE "LAST")
endif()

find_program(R_EXECUTABLE NAMES R R.exe)
#message(STATUS "R_EXECUTABLE=" ${R_EXECUTABLE})

#---searching R installation using R executable
if(R_EXECUTABLE)
  execute_process(COMMAND ${R_EXECUTABLE} RHOME
                  OUTPUT_VARIABLE R_ROOT_DIR
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  #message(STATUS "R_ROOT_DIR=" ${R_ROOT_DIR})
  
  find_path(R_INCLUDE_DIR R.h
            HINTS ${R_ROOT_DIR}
            PATHS /usr/local/lib /usr/local/lib64 /usr/share
            PATH_SUFFIXES include R/include
            DOC "Path to file R.h")
  #message(STATUS "R_INCLUDE_DIR=" ${R_INCLUDE_DIR})
  
  if (WIN32)
    find_file(R_LIBRARY R.dll
              HINTS ${R_ROOT_DIR}
              PATH_SUFFIXES lib bin/x64
              DOC "R library (R.dll).")
  else()
    find_library(R_LIBRARY R
                 HINTS ${R_ROOT_DIR}/lib
                 DOC "R library (example libR.a, libR.dylib, etc.).")
  endif()
  #message(STATUS "R_LIBRARY=" ${R_LIBRARY})
  
endif()

#---setting include dirs and libraries
set(R_LIBRARIES ${R_LIBRARY})
set(R_INCLUDE_DIRS ${R_INCLUDE_DIR})
foreach(_cpt ${R_FIND_COMPONENTS})
  execute_process(COMMAND echo "cat(find.package('${_cpt}'))"
                  COMMAND ${R_EXECUTABLE} --vanilla --slave
                  RESULT_VARIABLE _rc
                  ERROR_QUIET
                  OUTPUT_VARIABLE _cpt_path
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT _rc)
    set(R_${_cpt}_FOUND 1)
  endif()

  # TODO : This doesn't work under Windows, I guess.
  find_library(R_${_cpt}_LIBRARY
               lib${_cpt}.so lib${_cpt}.dylib
               HINTS ${_cpt_path}/lib)
  if(R_${_cpt}_LIBRARY)
    mark_as_advanced(R_${_cpt}_LIBRARY)
    list(APPEND R_LIBRARIES ${R_${_cpt}_LIBRARY})
  endif()

  find_path(R_${_cpt}_INCLUDE_DIR ${_cpt}.h HINTS  ${_cpt_path} PATH_SUFFIXES include R/include)
  if(R_${_cpt}_INCLUDE_DIR)
    mark_as_advanced(R_${_cpt}_INCLUDE_DIR)
    list(APPEND R_INCLUDE_DIRS ${R_${_cpt}_INCLUDE_DIR})
  endif()

endforeach()

# Handle the QUIETLY and REQUIRED arguments and set R_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(R HANDLE_COMPONENTS REQUIRED_VARS R_EXECUTABLE R_INCLUDE_DIR R_LIBRARY)
mark_as_advanced(R_FOUND R_EXECUTABLE R_INCLUDE_DIR R_LIBRARY)

# Find the R version and dispay message
execute_process(COMMAND ${R_EXECUTABLE} CMD BATCH --no-echo --no-timing ${CMAKE_CURRENT_SOURCE_DIR}/GetVersion.R ${CMAKE_CURRENT_BINARY_DIR}/GetVersion.Rout)
execute_process(COMMAND ${CMAKE_COMMAND} -E cat ${CMAKE_CURRENT_BINARY_DIR}/GetVersion.Rout
                ERROR_QUIET
                OUTPUT_VARIABLE R_VERSION
                OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "Found R: " ${R_LIBRARIES} " (found version \"" ${R_VERSION} "\")")



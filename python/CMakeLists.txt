if(NOT BUILD_PYTHON)
  return()
endif()

# Look for Python3
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

# Look for SWIG 4
find_package(SWIG 4.0 REQUIRED)
include(${SWIG_USE_FILE})

# Activate verbose mode
#set(CMAKE_VERBOSE_MAKEFILE ON)

# Display debug info for find_package
#set(CMAKE_FIND_DEBUG_MODE TRUE)

# Look for Boost dependencies
#find_package(Boost REQUIRED)

# Look for HDF5
# https://stackoverflow.com/questions/41529774/cmakelists-txt-for-compiling-hdf5
#find_package(HDF5 REQUIRED COMPONENTS C CXX)

# Look for pybind11
find_package(pybind11 QUIET) # try locally
if(NOT pybind11_FOUND) # else fetch it from github
  message(STATUS "Could NOT find pybind11 (fetching it from GitHub)")
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.6.1
  )
  FetchContent_MakeAvailable(pybind11)
endif()

# Create empty documentation
# TODO : put documentation.i out of the source
execute_process(
  COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
  COMMAND_ECHO STDOUT
)

  # Update documentation from gstlearn's XML
  # TODO : make this step optional (if XML doc exists only) and cross-platform
if (INCLUDE_DOXYGEN)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE}
                         ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxy2swig.py
                         ${PROJECT_BINARY_DIR}/doxygen/xml/index.xml
                         ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
    COMMAND_ECHO STDOUT
  )
endif()

# Keep sync with interface file content (*.i)
set(PYTHON_PACKAGE_NAME ${PROJECT_NAME})

# Generation folders
if (MSVC) # TODO: Prefer using GENERATOR_IS_MULTI_CONFIG
  set(PYTHON_PACKAGE_ROOT_FOLDER        ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>)
else()
  set(PYTHON_PACKAGE_ROOT_FOLDER        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
endif()
set(PYTHON_PACKAGE_DESTINATION_FOLDER ${PYTHON_PACKAGE_ROOT_FOLDER}/${PYTHON_PACKAGE_NAME})
# Location for SWIG target language code
set(CMAKE_SWIG_OUTDIR                 ${PYTHON_PACKAGE_DESTINATION_FOLDER})
# Location for SWIG wrapper library
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY    ${PYTHON_PACKAGE_DESTINATION_FOLDER}) # Use LIBRARY (gcc)!
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY    ${PYTHON_PACKAGE_DESTINATION_FOLDER}) # Use RUNTIME (msvc)!

# Create pygstlearn target
set(SRC pygstlearn.i)
set_source_files_properties(${SRC} PROPERTIES
  CPLUSPLUS ON
  # Swig wrapper will be regenerated if C++ headers are changed
  # https://stackoverflow.com/a/70375654/3952924
  #USE_SWIG_DEPENDENCIES TRUE   # Incompatible with CMAKE_SWIG_OUTDIR under cmake 3.22 (MSVC)
  # myfibo C++ headers will be found by %import directive
  USE_TARGET_INCLUDE_DIRECTORIES TRUE
  # Python package name will be the name of the installed folder (no module)
  SWIG_MODULE_NAME __init__
)

# Tell to CMake that ${SRC} depends on some other files
# TODO: This only works on GCC because MSVC always recompile swig (don't know why)
set_property(SOURCE ${SRC} PROPERTY
  DEPENDS
    ${PROJECT_SOURCE_DIR}/swig/swig_exp.i
    ${PROJECT_SOURCE_DIR}/swig/swig_inc.i
)

# Customize swig telling required Python 3+ version and keyword feature
# (WARNING: keep this BEFORE swig_add_library)
set(CMAKE_SWIG_FLAGS -py3 -keyword)
# Note: Cannot use -doxygen feature because all documentation is in C++ body files
# Note: Use -E option for seeing preprocessor output

# Add the swig target
swig_add_library(
  ${PYTHON_PACKAGE_NAME}
  TYPE SHARED
  LANGUAGE python
  SOURCES ${SRC}
)

# Link to static ${PROJECT_NAME} library
target_link_libraries(${PYTHON_PACKAGE_NAME} PUBLIC ${PROJECT_NAME}::static)

# Find appropriate Python.h
target_include_directories(${PYTHON_PACKAGE_NAME} PRIVATE ${Python3_INCLUDE_DIRS})

# Link to python's libraries (only for MSVC - Release)
if(MSVC)
  target_link_libraries(${PYTHON_PACKAGE_NAME} PRIVATE ${Python3_LIBRARIES})
endif()

# Link to pybind11 
target_link_libraries(${PYTHON_PACKAGE_NAME} PRIVATE pybind11::module) # TODO: use swig_link_library

# Rename the python module (without 'py' prefix)
#set_target_properties(${PYTHON_PACKAGE_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

# Generate setup.py automatically for each configuration
# First step: replace variables (@VAR@)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in 
  ${CMAKE_CURRENT_BINARY_DIR}/setup.py.in
  @ONLY
)
# Second step: replace generator expression ($<>)
file(GENERATE
  OUTPUT ${PYTHON_PACKAGE_ROOT_FOLDER}/setup.py
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py.in
)

# Generate README.md for each configuration
file(GENERATE
  OUTPUT ${PYTHON_PACKAGE_ROOT_FOLDER}/README.md
  INPUT  ${CMAKE_CURRENT_SOURCE_DIR}/README.md
)

# Copy each modules to destination folder
file(GLOB MODULES modules/*.py)
foreach(MODULE ${MODULES})
  file(COPY ${MODULE} DESTINATION ${PYTHON_PACKAGE_DESTINATION_FOLDER}/)
endforeach()

# Install python package (if SWIG target has run)
# TODO: Do also installation each time setup.py.in or README.md is modified
add_custom_command(TARGET ${PYTHON_PACKAGE_NAME}
  COMMAND ${Python3_EXECUTABLE} -m pip install ${PYTHON_PACKAGE_ROOT_FOLDER}
  COMMENT "Installing python package"
  VERBATIM
)
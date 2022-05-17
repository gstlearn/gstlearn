if(NOT BUILD_PYTHON)
  return()
endif()

# Python package name
set(PYTHON_PACKAGE_NAME ${PROJECT_NAME})
# Beware : the real package name is the PROJECT_NAME (for better import)

###################################
# DEPENDENCIES

# Look for Python3
#set(CMAKE_FIND_DEBUG_MODE TRUE)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)

if (APPLE)
  find_package(PythonLibs REQUIRED)
endif()

# Look for SWIG 4
find_package(SWIG 4.0 REQUIRED)
# TODO : If SWIG not found, fetch it from the web
include(${SWIG_USE_FILE})

# TODO : Look for matplotlib or fetch it from the web

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

######################################  
# DOCUMENTATION

# Ensure documentation.i is existing
# TODO : put documentation.i out of the source
add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
  COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
  COMMENT "Ensure that documentation.i exists")

# Generate python documentation (before wrapper)
# TODO : Do not regenerate python documentation if doxygen output is up to date
add_custom_target(python_doc
  COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxy2swig.py -q
    ${PROJECT_BINARY_DIR}/doxygen/xml/index.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
  COMMENT "Generate python documentation from doxygen output"
)

# Tell that we need to generate doxymentation before building python package
add_dependencies(python_doc doxygen)

######################################
# WRAPPER GENERATION & BUILD

# Python package interface file
set(SRC pygstlearn.i)
# No need to define %module in the interface file

# Generation folders
if (MSVC) # TODO: Prefer using GENERATOR_IS_MULTI_CONFIG ?
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

# SWIG interface source file properties
set_source_files_properties(${SRC} PROPERTIES
  CPLUSPLUS ON
  # Swig wrapper will be regenerated if C++ headers are changed
  # https://stackoverflow.com/a/70375654/3952924
  USE_SWIG_DEPENDENCIES TRUE
  # gstlearn C++ headers will be found by %import directive
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
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
)

# SWIG flags
set(CMAKE_SWIG_FLAGS -py3 -keyword)
# Note: Cannot use -doxygen feature because all documentation is in C++ body files
# Note: Use -E option for seeing SWIG preprocessor output

# SWIG target (generate wrapper & python script then build [so,pyd] library)
swig_add_library(
  python_build
  LANGUAGE   python
  TYPE       SHARED
  SOURCES    ${SRC}
)

# Link to static library
target_link_libraries(python_build PUBLIC ${PROJECT_NAME}::static)

# Transmit the fact that we are linking to static ${PROJECT_NAME} library
if (MSVC)
  # Wrapper source is too big for MSVC !
  set(COMP_FLAGS "/bigobj -D${PROJECT_NAME_UP}_STATIC_DEFINE")
else()
  set(COMP_FLAGS "-D${PROJECT_NAME_UP}_STATIC_DEFINE")
endif()

# Set some properties on SWIG target
set_target_properties(python_build PROPERTIES
  # Do not build python package by default
  EXCLUDE_FROM_ALL ON
  # Modify the generated library name
  OUTPUT_NAME ${PROJECT_NAME}
  # Compiler flags for the SWIG library
  COMPILE_FLAGS ${COMP_FLAGS}
)

# Find appropriate Python.h
target_include_directories(python_build PRIVATE ${Python3_INCLUDE_DIRS})

# Link to python's libraries (only for MSVC and MacOS - Release)
if(MSVC)
  target_link_libraries(python_build PRIVATE ${Python3_LIBRARIES})
endif()

if(APPLE)
  target_link_libraries(python_build PRIVATE ${PYTHON_LIBRARIES})
endif()


# Link to pybind11 
target_link_libraries(python_build PRIVATE pybind11::module) # TODO: use swig_link_library

########################################
# PACKAGING (target independent)

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

# Copy each modules to package *destination* folder for each configuration
file(GLOB MODULES ${CMAKE_CURRENT_SOURCE_DIR}/modules/*.py)
foreach(MODULE ${MODULES})
  cmake_path(GET MODULE FILENAME MODULE_NAME)
  set(MODULE_NAME ${PYTHON_PACKAGE_NAME}/${MODULE_NAME})
  file(GENERATE
    OUTPUT ${PYTHON_PACKAGE_ROOT_FOLDER}/${MODULE_NAME}
    INPUT ${MODULE}
  )
endforeach()

######################################
# INSTALLATION

# Add a custom target for python package installation
# TODO: Do also installation each time setup.py.in or README.md is modified
add_custom_target(python_install
  COMMAND ${Python3_EXECUTABLE} -m pip install ${PYTHON_PACKAGE_ROOT_FOLDER}
  COMMENT "Installing python package"
  VERBATIM
)

# Tell to CMake that python package must be built before installation [python_doc is optional]
add_dependencies(python_install python_build)

######################################
# UPLOAD

# Warning : Need to customize your TestPyPi token :
# https://pypi.org/help/#apitoken
 
# Add a custom target for python package distribution uploading to PyPi
# TODO: Do also installation each time setup.py.in or README.md is modified
if (MSVC) # Wheel binary package
  add_custom_target(python_upload
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${PYTHON_PACKAGE_ROOT_FOLDER}/dist
    COMMAND ${Python3_EXECUTABLE} setup.py bdist_wheel
#    COMMAND ${Python3_EXECUTABLE} -m twine upload --repository testpypi dist/*.whl
    COMMENT "Building and uploading python package distribution"
    WORKING_DIRECTORY ${PYTHON_PACKAGE_ROOT_FOLDER}
    VERBATIM
  )
else() # Source package # TODO : To be confirmed (because gstlearn.a)
  add_custom_target(python_upload
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${PYTHON_PACKAGE_ROOT_FOLDER}/dist
    COMMAND ${Python3_EXECUTABLE} setup.py sdist
#    COMMAND ${Python3_EXECUTABLE} -m twine upload --repository testpypi dist/*.tar.gz
    COMMENT "Building and uploading python package distribution"
    WORKING_DIRECTORY ${PYTHON_PACKAGE_ROOT_FOLDER}
    VERBATIM
  )
endif()

# Tell to CMake that python package must be built before uploading [python_doc is optional]
add_dependencies(python_upload python_build)

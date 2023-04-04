if(NOT BUILD_PYTHON)
  return()
endif()

# Python package name
set(PYTHON_PACKAGE_NAME ${PROJECT_NAME})
# Beware : the real package name is the PROJECT_NAME (for better import)

###################################
# DEPENDENCIES

# Look for Python3
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# Look for SWIG 4 [No need for customized swig]
find_package(SWIG 4.0 REQUIRED)
message(STATUS "Found SWIG: " ${SWIG_EXECUTABLE} " (found version \"" ${SWIG_VERSION} "\")")
include(${SWIG_USE_FILE})

# Look for pybind11
find_package(pybind11 QUIET) # try locally
if(NOT pybind11_FOUND) # else fetch it from github (keep this as long as runner have no pybind11-global)
  message(STATUS "Could NOT find pybind11 (fetching it from GitHub)")
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    #GIT_TAG v2.6.1 # Comment: Use last pybind11 version available
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

# Doxygen is optional
if (BUILD_DOXYGEN)
  add_dependencies(python_doc doxygen)
endif()

######################################
# WRAPPER GENERATION & BUILD

# Python package interface file
set(SRC pygstlearn.i)

# Generation folders
if (IS_MULTI_CONFIG)
  set(PYTHON_PACKAGE_ROOT_FOLDER      ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>)
else()
  set(PYTHON_PACKAGE_ROOT_FOLDER      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
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
  # https://stackoverflow.com/a/70375654/3952924 (CMake 3.20 minimum)
  USE_SWIG_DEPENDENCIES TRUE
  # C++ headers will be found by %import swig directive
  USE_TARGET_INCLUDE_DIRECTORIES TRUE
  # Python package name will be the name of the installed folder (no module)
  SWIG_MODULE_NAME __init__
)

# Tell to CMake that ${SRC} depends on some other files
set_property(SOURCE ${SRC} PROPERTY
  DEPENDS
    ${PROJECT_SOURCE_DIR}/swig/swig_exp.i
    ${PROJECT_SOURCE_DIR}/swig/swig_inc.i
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
)

# SWIG flags
set(CMAKE_SWIG_FLAGS -keyword) # Use -keyword option to activate kwargs
# Note: Cannot use -doxygen feature because all documentation is in C++ body files
# Note: Use -E option for seeing SWIG preprocessor output

# SWIG target (generate wrapper & python script then build [so,pyd] library)
swig_add_library(
  python_build
  LANGUAGE   python
  TYPE       SHARED
  SOURCES    ${SRC}
)

if (APPLE)
  # To prevent 'import' from crashing under Mac (swig output suffix is .dynlib by default) 
  set_target_properties(python_build PROPERTIES SUFFIX ".so")
endif()

# Link to static library and NumPy
target_link_libraries(python_build PUBLIC ${PROJECT_NAME}::static Python3::NumPy)

# Transmit the fact that we are linking to static ${PROJECT_NAME} library
set(COMP_FLAGS "-D${PROJECT_NAME_UP}_STATIC_DEFINE")
if (MSVC)
  set(COMP_FLAGS "${COMP_FLAGS} /wd4706 /wd4456 /wd4127 /wd4100") # SWIG generated code warnings
endif()
if (MSVC)
  # Wrapper source is too big for MSVC !
  set(COMP_FLAGS "${COMP_FLAGS} /bigobj")
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

if (MSVC OR APPLE)
  # Link to python's libraries (needed for MSVC and MacOS)
  # Prefer using Python3::Module (which works for MacOS) in place of ${Python3_LIBRARIES}
  # https://gitlab.kitware.com/cmake/cmake/-/issues/17664
  target_link_libraries(python_build PRIVATE Python3::Module)
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

# Generate version.py automatically for each configuration
# First step: replace variables (@VAR@)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.py.in 
  ${CMAKE_CURRENT_BINARY_DIR}/version.py.in
  @ONLY
)

# Second step: replace generator expression ($<>)
file(GENERATE
  OUTPUT ${PYTHON_PACKAGE_ROOT_FOLDER}/${PYTHON_PACKAGE_NAME}/version.py
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/version.py.in
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

# Warning: Need to customize your TestPyPi token:
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

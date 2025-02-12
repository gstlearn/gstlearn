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
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module NumPy)
string(REGEX MATCH "^([0-9]+)\\.([0-9]+)\\.([0-9]+)"
  Python3_NumPy_VERSION_MATCH ${Python3_NumPy_VERSION})
set(Python3_NumPy_VERSION_MAJOR ${CMAKE_MATCH_1})
math(EXPR Python3_NumPy_NEXT_MAJOR "${Python3_NumPy_VERSION_MAJOR} + 1")

# Look for SWIG 4 [No need for customized swig]
find_package(SWIG 4.0 REQUIRED)
include(${SWIG_USE_FILE})

######################################  
# DOCUMENTATION

# Ensure documentation.i is existing
# TODO : put documentation.i out of the source (be careful of pygstlearn.i)
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

####CUSTOM WRAPPER
add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/python/generated_python.i
  COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/python/generated_python.i
  COMMENT "Ensure that generated_python.i exists")
  
add_custom_target(python_wrap_macro
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../swig/macrogenerate.py
  "python"
  COMMENT "Generate macro wrapper package"
  )


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

# Reuse shared library objects
target_sources(
  python_build
  PRIVATE
  $<TARGET_OBJECTS:${PROJECT_NAME}::shared>
)
target_include_directories(python_build
  PRIVATE
  $<TARGET_PROPERTY:${PROJECT_NAME}::shared,INCLUDE_DIRECTORIES>
)
target_link_libraries(python_build
  PRIVATE
  $<TARGET_PROPERTY:${PROJECT_NAME}::shared,LINK_LIBRARIES>
)

# Link to NumPy
target_link_libraries(python_build PUBLIC Python3::NumPy)

# Transmit the fact that we are linking to static ${PROJECT_NAME} library
set(COMP_FLAGS "-D${PROJECT_NAME_UP}_STATIC_DEFINE")
if (MSVC)
  set(COMP_FLAGS "${COMP_FLAGS} /wd4706 /wd4456 /wd4127 /wd4100") # SWIG generated code warnings
else ()
  set(COMP_FLAGS "${COMP_FLAGS} -Wno-unused-parameter -Wno-free-nonheap-object") # Prevent some SWIG C++/Python warnings
endif()
if (MSVC)
  # Wrapper source is too big for MSVC !
  set(COMP_FLAGS "${COMP_FLAGS} /bigobj")
endif()

# Set some properties on SWIG target
set_target_properties(python_build PROPERTIES
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

if(USE_BOOST_SPAN)
  target_compile_definitions(python_build PUBLIC USE_BOOST_SPAN)
endif()

add_dependencies(python_build python_wrap_macro)
add_dependencies(python_build generate_output)

########################################
# PACKAGING (target independent)

# Generate setup.py automatically for each configuration
# First step: replace variables (@VAR@)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/setup.py.in
  @ONLY
)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml.in
  ${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml.in
  @ONLY
)

# Second step: replace generator expression ($<>)
file(GENERATE
  OUTPUT ${PYTHON_PACKAGE_ROOT_FOLDER}/setup.py
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py.in
)
file(GENERATE
  OUTPUT ${PYTHON_PACKAGE_ROOT_FOLDER}/pyproject.toml
  INPUT ${CMAKE_CURRENT_BINARY_DIR}/pyproject.toml.in
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
  COMMAND ${Python3_EXECUTABLE} -m pip install "${PYTHON_PACKAGE_ROOT_FOLDER}[plot]"
  COMMENT "Installing python package"
  )
 
# Tell to CMake that python package must be built before installation [python_doc is optional]

add_dependencies(python_install python_build)

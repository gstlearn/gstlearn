# R package name
set(R_PACKAGE_NAME ${PROJECT_NAME})
# Beware : the real package name is the PROJECT_NAME (for better import)

###################################
# DEPENDENCIES

# Look for R
include(FindR.cmake)

# Look for SWIG 4
find_package(SWIG 4.0 REQUIRED)
# TODO : If SWIG not found, fetch it from the web
include(${SWIG_USE_FILE})

######################################  
# DOCUMENTATION
#
# Ensure documentation.i is existing
# TODO : put documentation.i out of the source
#add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
#  COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
#  COMMENT "Ensure that documentation.i exists")
#
# Generate python documentation (before wrapper)
# TODO : Do not regenerate python documentation if doxygen output is up to date
#add_custom_target(python_doc
#  COMMAND ${Python3_EXECUTABLE}
#    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxy2swig.py -q
#    ${PROJECT_BINARY_DIR}/doxygen/xml/index.xml
#    ${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
#  COMMENT "Generate python documentation from doxygen output"
#)
#
# Tell that we need to generate doxymentation before building python package
#add_dependencies(python_doc doxygen)

######################################
# WRAPPER GENERATION & BUILD

# R package target
set(RTGT r_build)

# R package interface file
set(SRC rgstlearn.i)

# Generation folders
if (MSVC) # TODO: Prefer using GENERATOR_IS_MULTI_CONFIG
  set(R_PACKAGE_ROOT_FOLDER        ${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>)
else()
  set(R_PACKAGE_ROOT_FOLDER        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
endif()
set(R_PACKAGE_DESTINATION_FOLDER ${R_PACKAGE_ROOT_FOLDER}/${R_PACKAGE_NAME})
# Location for SWIG target language code
set(CMAKE_SWIG_OUTDIR                 ${R_PACKAGE_DESTINATION_FOLDER})
# Location for SWIG wrapper library
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY    ${R_PACKAGE_DESTINATION_FOLDER}) # Use LIBRARY (gcc)!
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY    ${R_PACKAGE_DESTINATION_FOLDER}) # Use RUNTIME (msvc)!

# SWIG interface source file properties
set_source_files_properties(${SRC} PROPERTIES
  CPLUSPLUS ON
  # Swig wrapper will be regenerated if C++ headers are changed
  # https://stackoverflow.com/a/70375654/3952924
  USE_SWIG_DEPENDENCIES TRUE
  # gstlearn C++ headers will be found by %import directive
  USE_TARGET_INCLUDE_DIRECTORIES TRUE
)

# Tell to CMake that ${SRC} depends on some other files
# TODO: This only works on GCC because MSVC always recompile swig (don't know why)
set_property(SOURCE ${SRC} PROPERTY
  DEPENDS
    ${PROJECT_SOURCE_DIR}/swig/swig_exp.i
    ${PROJECT_SOURCE_DIR}/swig/swig_inc.i
    #${CMAKE_CURRENT_SOURCE_DIR}/doc/documentation.i
)

# SWIG flags
#set(CMAKE_SWIG_FLAGS)
# Note: Cannot use -doxygen feature because all documentation is in C++ body files
# Note: Use -E option for seeing SWIG preprocessor output

# SWIG target (generate wrapper & R script then build [so,dll] library)
swig_add_library(
  ${RTGT}
  LANGUAGE   r
  TYPE       SHARED
  SOURCES    ${SRC}
)

# Link to static library
target_link_libraries(${RTGT} PUBLIC ${PROJECT_NAME}::static)

# Transmit the fact that we are linking to static ${PROJECT_NAME} library
if (MSVC)
  # Wrapper source is too big for MSVC !
  set(COMP_FLAGS "/bigobj -D${PROJECT_NAME_UP}_STATIC_DEFINE")
else()
  set(COMP_FLAGS "-D${PROJECT_NAME_UP}_STATIC_DEFINE")
endif()
# Set some properties on SWIG target
set_target_properties(${RTGT} PROPERTIES
  # Do not build R package by default
  EXCLUDE_FROM_ALL ON
  # Modify the generated library name
  OUTPUT_NAME ${PROJECT_NAME}
  # Compiler flags for the SWIG library
  COMPILE_FLAGS ${COMP_FLAGS}
)

# Find appropriate R.h
target_include_directories(${RTGT} PRIVATE ${R_INCLUDE_DIRS})

# Link to python's libraries (only for MSVC - Release)
if(MSVC)
  target_link_libraries(${RTGT} PRIVATE ${R_LIBRARIES})
endif()

if(APPLE)
  target_link_libraries(${RTGT} PRIVATE ${R_LIBRARIES})
endif()



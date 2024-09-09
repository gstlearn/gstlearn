# Make Release version the default (only for single configuration generators)
# TODO : Differentiate build directories for Debug and Release
if(NOT IS_MULTI_CONFIG)
  if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
  endif()
  # Show current configuration
  message(STATUS "BUILD_TYPE=" ${CMAKE_BUILD_TYPE})
endif()

# Add c++11 support whatever the compiler
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Warning fiesta!
# https://cmake.org/cmake/help/latest/command/add_compile_options.html
if (MSVC)
  # Warning level 4 (4 = maximum, 0 = none)
  add_compile_options(/bigobj /W4 /wd4251 /wd4244) # Except those two warnings
  # Silence MSVC warnings about unsafe C standard library functions
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
  # Lots of warnings (-Wall = add some warnings, -Wextra = add a ton of warnings)
  add_compile_options(
    -Wall
    -Wextra
    -Wno-deprecated-copy
    -Wtype-limits
    -Wnon-virtual-dtor
    -Wvla
    -Wundef
  )
  if (APPLE)
    add_compile_options(-Wno-absolute-value -Wno-inconsistent-missing-override)
  endif()
endif()

# C++ header location (keep the trailing '/')
set(INCLUDES 
    ${PROJECT_SOURCE_DIR}/include/)
# C++ source path (prevent using GLOB)
include(src/all_sources.cmake)
set(SOURCES)
foreach(CPP ${SRC})
  set(SOURCES ${SOURCES} ${PROJECT_SOURCE_DIR}/src/${CPP})
endforeach(CPP ${SRC})

# target_include_directories() below is enough to get the code to compile but if
# includes are not explicitly added to SOURCES then Visual Studio doesn't see them
# as part of the project itself. They are GLOB'ed even though it's fragile for
# simplicity (it's just to get the project tree right).
if (CMAKE_GENERATOR MATCHES "Visual Studio")
  file(GLOB_RECURSE INCS RELATIVE ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/include/*.hpp ${PROJECT_SOURCE_DIR}/include/*.h)
  foreach(HPP ${INCS})
    set(SOURCES ${SOURCES} ${PROJECT_SOURCE_DIR}/include/${HPP})
  endforeach(HPP ${INCS})
endif()

# Generation folder (into Release or Debug)
if (NOT IS_MULTI_CONFIG AND NOT WIN32)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE} OUTPUT_VARIABLE CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE} OUTPUT_VARIABLE CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE} OUTPUT_VARIABLE CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
endif()

# Debug find package instruction
#set(CMAKE_FIND_DEBUG_MODE TRUE)

# Look for Boost
#set(Boost_DEBUG 1)
find_package(Boost REQUIRED)
# TODO : If Boost not found, fetch it from the web ?

# Look for OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  add_definitions(-DOPENMP)
endif()

# Look for Eigen
find_package(Eigen3 REQUIRED) 
if(Eigen3_FOUND) 
    message(STATUS "Found Eigen3 version ${Eigen3_VERSION} in ${Eigen3_DIR}")
endif()

# Look for HDF5
if (USE_HDF5)
  # Use static library for HDF5 under Windows (no more issue with DLL location)
  if (WIN32)
    set(HDF5_USE_STATIC_LIBRARIES ON)
  endif()
  
  # Look for HDF5
  # https://stackoverflow.com/questions/41529774/cmakelists-txt-for-compiling-hdf5
  find_package(HDF5 REQUIRED COMPONENTS CXX)
  # TODO : If HDF5 not found, fetch it from the web ?
endif()

# Shared and Static libraries
add_library(shared                  SHARED ${SOURCES})
add_library(static EXCLUDE_FROM_ALL STATIC ${SOURCES})
set(FLAVORS shared static)

############################## Loop on flavors: shared and static
foreach(FLAVOR ${FLAVORS})
  # Convert flavor to uppercase
  string(TOUPPER ${FLAVOR} FLAVOR_UP)

  # Alias target for a better name
  add_library(${PROJECT_NAME}::${FLAVOR} ALIAS ${FLAVOR})

  # Include directories
  # PUBLIC is mandatory for tests and packages (no need to install)
  target_include_directories(${FLAVOR} PUBLIC
    # Add includes path for compiling the library
    $<BUILD_INTERFACE: ${INCLUDES}>
    # Add binary directory to find generated version.h and export.hpp
    $<BUILD_INTERFACE: ${PROJECT_BINARY_DIR}>
  )

  # Set some target properties
  set_target_properties(${FLAVOR} PROPERTIES
    # Hide all symbols by default (impose same behavior between Linux and Windows)
    C_VISIBILITY_PRESET hidden
    CXX_VISIBILITY_PRESET hidden
    # Any client who links the library needs -fPIC (static or shared)
    POSITION_INDEPENDENT_CODE 1
  )

  # Rename the output library name
  set_target_properties(${FLAVOR} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  # append a 'd' to the output file name of the debug build targets
  set_target_properties(${FLAVOR} PROPERTIES DEBUG_POSTFIX "d")
  
  # Set library version
  set_target_properties(${FLAVOR} PROPERTIES VERSION ${PROJECT_VERSION})
  
  # Enable OpenMP
  target_link_libraries(${FLAVOR} PRIVATE OpenMP::OpenMP_CXX)

  # Link to csparse and gmtsph
  target_link_libraries(${FLAVOR} PRIVATE csparse gmtsph)
    
  # Link to Eigen
  target_link_libraries(${FLAVOR} PUBLIC Eigen3::Eigen)

  # Link to Boost (use headers)
  # Target for header-only dependencies. (Boost include directory)
  target_link_libraries(${FLAVOR} PRIVATE Boost::boost)
  
  # Link to HDF5
  if (USE_HDF5)
    # Define _USE_HDF5 macro
    target_compile_definitions(${FLAVOR} PUBLIC _USE_HDF5) 
    target_link_libraries(${FLAVOR} PUBLIC hdf5::hdf5_cpp)
  endif()
  
  # Exclude [L]GPL features from Eigen
  #target_compile_definitions(${FLAVOR} PUBLIC EIGEN_MPL2_ONLY) 

  # Link to specific libraries (only for Microsoft Visual Studio)
  if (MSVC)
    target_link_libraries(${FLAVOR} PUBLIC iphlpapi rpcrt4)
  endif()
  if (MINGW)
    target_link_libraries(${FLAVOR} PUBLIC -liphlpapi -lrpcrt4)
  endif()

  # Build a cmake file to be imported by library users
  export(TARGETS ${FLAVOR}
         NAMESPACE ${PROJECT_NAME}::
         FILE ${GSTLEARN_CMAKE_FILE}
         APPEND)
         
endforeach(FLAVOR ${FLAVORS})
############################## End loop on flavors


###################### Shared library specific options

# Generate export header
include(GenerateExportHeader)
set(DISABLE_EXPORT_IF_SWIG "
#ifdef SWIG
#    undef ${PROJECT_NAME_UP}_EXPORT
#    undef ${PROJECT_NAME_UP}_NO_EXPORT
#    undef ${PROJECT_NAME_UP}_DEPRECATED
#    undef ${PROJECT_NAME_UP}_DEPRECATED_EXPORT
#    undef ${PROJECT_NAME_UP}_DEPRECATED_NO_EXPORT
#    define ${PROJECT_NAME_UP}_EXPORT
#    define ${PROJECT_NAME_UP}_NO_EXPORT
#    define ${PROJECT_NAME_UP}_DEPRECATED
#    define ${PROJECT_NAME_UP}_DEPRECATED_EXPORT
#    define ${PROJECT_NAME_UP}_DEPRECATED_NO_EXPORT
#endif
#ifdef ${PROJECT_NAME_UP}_STATIC_DEFINE
#    define ${PROJECT_NAME_UP}_TEMPLATE_EXPORT
#else
#    ifdef shared_EXPORTS
#        define ${PROJECT_NAME_UP}_TEMPLATE_EXPORT
#    else
#        define ${PROJECT_NAME_UP}_TEMPLATE_EXPORT extern
#    endif
#endif
")
generate_export_header(shared
  BASE_NAME ${PROJECT_NAME}
  EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/${PROJECT_NAME}_export.hpp
  CUSTOM_CONTENT_FROM_VARIABLE DISABLE_EXPORT_IF_SWIG
)

# Set the so version to project major version
set_target_properties(shared PROPERTIES
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

###################### Static library specific options

# Prevent from using _declspec when static
set_target_properties(static PROPERTIES
  COMPILE_FLAGS -D${PROJECT_NAME_UP}_STATIC_DEFINE
)

# we need a specific name for the static library otherwise Ninja on
# Windows not happy...
if (WIN32 AND CMAKE_GENERATOR MATCHES "Ninja")
  set_target_properties(static PROPERTIES OUTPUT_NAME ${PROJECT_NAME}_static)
endif()

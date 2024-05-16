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
else()
  # Lots of warnings (-Wall = add some warnings, -Wextra = add a ton of warnings)
  add_compile_options(-Wall -Wextra -Wno-deprecated-copy -Wno-unused-parameter)
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

# Generation folder (into Release or Debug)
if (NOT IS_MULTI_CONFIG)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE} OUTPUT_VARIABLE CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE} OUTPUT_VARIABLE CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE} OUTPUT_VARIABLE CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
endif()

# Change the name of the output file (to distinguish lib files under Windows)
if (WIN32)
  set(CMAKE_STATIC_LIBRARY_PREFIX "lib")
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
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#  if(${APPLE})
#    include_directories(${OpenMP_C_INCLUDE_DIR})
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lomp") # clang++: warning: -lomp: 'linker' input unused
#  endif()
endif()

# Look for Eigen
find_package(Eigen3 REQUIRED) 
if(EIGEN3_FOUND)
  file(READ "${EIGEN3_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" eigen_macros_h)
  string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" eigen_world "${eigen_macros_h}")
  set(EIGEN_VERSION_WORLD "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" eigen_major "${eigen_macros_h}")
  set(EIGEN_VERSION_MAJOR "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" eigen_minor "${eigen_macros_h}")
  set(EIGEN_VERSION_MINOR "${CMAKE_MATCH_1}")
  set(EIGEN_VERSION ${EIGEN_VERSION_WORLD}.${EIGEN_VERSION_MAJOR}.${EIGEN_VERSION_MINOR})
  message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR} (found version ${EIGEN_VERSION})")
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
    # Add Eigen include directories
    $<BUILD_INTERFACE: ${EIGEN3_INCLUDE_DIR}>
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
  
  # Link to csparse and gmtsph
  target_link_libraries(${FLAVOR} PRIVATE csparse gmtsph)
    
  # Link to Boost (use headers)
  # Target for header-only dependencies. (Boost include directory)
  target_link_libraries(${FLAVOR} PRIVATE Boost::boost)
  
  # Link to HDF5
  if (USE_HDF5)
    # Define _USE_HDF5 macro
    target_compile_definitions(${FLAVOR} PUBLIC _USE_HDF5) 
    
    # CMake>=3.19 introduces hdf5 targets that could be used the same way as boost targets
    target_include_directories(${FLAVOR} PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(${FLAVOR} PUBLIC ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
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

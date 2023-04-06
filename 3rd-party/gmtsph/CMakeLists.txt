# Add c++11 support whatever the compiler
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# C++ header location (keep the trailing '/')
set(INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/)
# C++ source path
set(SOURCES gmtsph.cpp)

# Change the name of the output file (to distinguish lib files under Windows)
#if (WIN32)
#  set(CMAKE_STATIC_LIBRARY_PREFIX "lib")
#endif()

# Static library
add_library(gmtsph STATIC ${SOURCES})

# Include directories
target_include_directories(gmtsph PUBLIC
  # Add includes path for compiling the library
  $<BUILD_INTERFACE: ${INCLUDES}>
)

# Set some target properties
set_target_properties(gmtsph PROPERTIES
  # Any client who links the library needs -fPIC
  POSITION_INDEPENDENT_CODE 1
)

# Export target
export(TARGETS gmtsph
       NAMESPACE ${PROJECT_NAME}::
       FILE ${GSTLEARN_CMAKE_FILE}
       APPEND)
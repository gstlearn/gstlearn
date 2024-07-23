# Add c++11 support whatever the compiler
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# C++ header location (keep the trailing '/')
set(INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/)
# C++ source path
set(SOURCES csparse.cpp)

# Static library
add_library(csparse STATIC ${SOURCES})

# Include directories
target_include_directories(csparse PUBLIC
  # Add includes path for compiling the library
  $<BUILD_INTERFACE: ${INCLUDES}>
)

# Set some target properties
set_target_properties(csparse PROPERTIES
  # Any client who links the library needs -fPIC
  POSITION_INDEPENDENT_CODE 1
)

# Export target
export(TARGETS csparse
       NAMESPACE ${PROJECT_NAME}::
       FILE ${GSTLEARN_CMAKE_FILE}
       APPEND)
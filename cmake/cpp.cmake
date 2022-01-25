# C++ code location
include(src/all_sources.cmake)
set(SOURCES)
foreach(CPP ${SRC})
  set(SOURCES ${SOURCES} ${PROJECT_SOURCE_DIR}/src/${CPP})
endforeach(CPP ${SRC})

# Generation folder
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE})

# Impose 'd' suffix in debug (global property)
set(CMAKE_DEBUG_POSTFIX d)

# Look for Boost
find_package(Boost REQUIRED)

# Look for HDF5
# https://stackoverflow.com/questions/41529774/cmakelists-txt-for-compiling-hdf5
find_package(HDF5 REQUIRED COMPONENTS C CXX)

# Shared library if needed
set(FLAVORS)
if (BUILD_SHARED)
  add_library(shared SHARED ${SOURCES})
  set(FLAVORS ${FLAVORS} shared)
endif()
  
# Static library if needed
if (BUILD_STATIC)
  add_library(static STATIC ${SOURCES})
  set(FLAVORS ${FLAVORS} static)
endif()

############################## Loop on flavor: shared and static
foreach(FLAVOR ${FLAVORS})
  # Convert flavor to uppercase
  string(TOUPPER ${FLAVOR} FLAVOR_UP)
  
  # Alias target for a better name
  add_library(${PROJECT_NAME}::${FLAVOR} ALIAS ${FLAVOR})
    
  # Include directories
  target_include_directories(${FLAVOR} PUBLIC
    # Includes for compiling the library
    $<BUILD_INTERFACE: ${PROJECT_SOURCE_DIR}/include>
    # Add binary directory to find generated version.h
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  )

  # Set some target properties
  set_target_properties(${FLAVOR} PROPERTIES
    # Symbol for export.hpp header (do not use GenerateExportHeader)
    COMPILE_FLAGS "-D${PROJECT_NAME_UP}_BUILD_${FLAVOR_UP}"
    # Hide all symbols by default (impose same behavior between Linux and Windows)
    C_VISIBILITY_PRESET hidden
    CXX_VISIBILITY_PRESET hidden
    # Any client who links the library needs -fPIC (static or shared)
    POSITION_INDEPENDENT_CODE 1
  )
  
  # Rename the output library name
  set_target_properties(${FLAVOR} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  
  # Set library version
  set_target_properties(${FLAVOR} PROPERTIES VERSION ${PROJECT_VERSION})
  
  # Used by delaunay (On Windows, prevent the include sys/time.h (-DNO_TIMER))
  target_compile_definitions(${FLAVOR} PRIVATE NO_TIMER) 
  
  # Link to boost
  # Target for header-only dependencies. (Boost include directory)
  # It should be PRIVATE if no headers of the gstlearn include boost files
  target_link_libraries(${FLAVOR} PUBLIC Boost::boost)
  
  # Link to HDF5
  # CMake>=3.19 introduces hdf5 targets that could be used the same way as boost targets
  target_include_directories(${FLAVOR} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${FLAVOR} PUBLIC ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

  # Link to specific libraries (only for Windows)
  if (WIN32)
    target_link_libraries(${FLAVOR} PUBLIC iphlpapi rpcrt4)
  endif()

endforeach(FLAVOR ${FLAVORS})
############################## End loop on flavor

# Shared library specific options
if (BUILD_SHARED)
  # Set the so version to project major version
  set_target_properties(shared PROPERTIES
    SOVERSION ${PROJECT_VERSION_MAJOR}
  )
endif()

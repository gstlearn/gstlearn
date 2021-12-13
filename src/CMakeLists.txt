# Include all sources
include(all_sources.cmake)

# Look for Boost
find_package(Boost REQUIRED)

# Look for HDF5
# https://stackoverflow.com/questions/41529774/cmakelists-txt-for-compiling-hdf5
find_package(HDF5 REQUIRED COMPONENTS C CXX)

# Define the targets
add_library(shared SHARED ${SOURCES})
add_library(static STATIC ${SOURCES})

############################## Loop on flavor: shared and static
set(FLAVORS shared static)
foreach(FLAVOR ${FLAVORS})
  
  # On Windows, prevent the include sys/time.h (-DNO_TIMER)
  target_compile_definitions(${FLAVOR} PRIVATE NO_TIMER)
  
  # Add gstlearn includes (add swig include directories (installed) for swig module who wants to link with gstlearn)
  target_include_directories(${FLAVOR} PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gstlearn>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/swig>"
  )
  
  # Add binary directory to find generated version.h and gstlearn_export.h
  target_include_directories(${FLAVOR} PUBLIC "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>")
  
  # Rename the output library name
  set_target_properties(${FLAVOR} PROPERTIES OUTPUT_NAME gstlearn) # Will produce libgstlearn.so (on Linux)
  
  # Set library version
  set_target_properties(${FLAVOR} PROPERTIES VERSION ${PROJECT_VERSION})
  
  # Anyone who wants to build a binary which links the gstlearn needs -fPIC
  set_property(TARGET ${FLAVOR} PROPERTY POSITION_INDEPENDENT_CODE 1)
  
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
    target_link_libraries(${FLAVOR} PUBLIC -liphlpapi -lrpcrt4)
  endif()

endforeach(FLAVOR ${FLAVORS})

############################## End loop on flavor

# Set the so version to project major version
set_target_properties(shared PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})

# Generate export header
include(GenerateExportHeader)
set(DISABLE_EXPORT_IF_SWIG "
#ifdef SWIG
#undef GSTLEARN_EXPORT
#undef GSTLEARN_NO_EXPORT
#define GSTLEARN_EXPORT
#define GSTLEARN_NO_EXPORT
#endif
")
generate_export_header(shared BASE_NAME gstlearn
                       EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/gstlearn_export.hpp
                       CUSTOM_CONTENT_FROM_VARIABLE DISABLE_EXPORT_IF_SWIG
)
set_target_properties(static PROPERTIES COMPILE_FLAGS -DGSTLEARN_STATIC_DEFINE)

# Install the gstlearn_export file
install(FILES ${PROJECT_BINARY_DIR}/gstlearn_export.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gstlearn
)

# Install the version file
install(FILES ${PROJECT_BINARY_DIR}/version.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gstlearn
)

# Install dynamic library in DESTINATION/lib folder
install(TARGETS shared EXPORT corelibs
        LIBRARY DESTINATION lib
        OPTIONAL
)
# Install static library in DESTINATION/lib folder
install(TARGETS static EXPORT corelibs
        ARCHIVE DESTINATION lib
        OPTIONAL
)

# Install the includes  
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/                 # Includes from source folder
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gstlearn       # in DESTINATION/include/gstlearn
)
# Install swig interface files
install(DIRECTORY ${CMAKE_SOURCE_DIR}/swig/                    # Swig files from source folder
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/swig # in DESTINATION/share/gstlearn/swig
)

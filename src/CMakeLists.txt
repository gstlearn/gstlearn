include(all_sources.cmake)

# Look for Boost
#option(Boost_DEBUG "Set to ON to enable debug output from FindBoost." ON)
#option(Boost_DETAILED_FAILURE_MSG "Set to ON to get detailed information" ON)
find_package(Boost REQUIRED)

# Look for HDF5
# https://stackoverflow.com/questions/41529774/cmakelists-txt-for-compiling-hdf5
find_package(HDF5 REQUIRED COMPONENTS C CXX)

############################## TODO : Start loop on flavor: shared and static

# Define targets 
add_library(shared SHARED ${SOURCES})
add_library(static STATIC ${SOURCES})

# On Windows, prevent the include sys/time.h (-DNO_TIMER)
target_compile_definitions(shared PRIVATE NO_TIMER)
target_compile_definitions(static PRIVATE NO_TIMER) 

# Add gstlearn includes (add swig include directories (installed) for swig module who wants to link with gstlearn)
target_include_directories(shared PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gstlearn>"   # TODO : still needed ?
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/swig>"
)
target_include_directories(static PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gstlearn>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/swig>"
)


# Add HDF5 include directory
target_include_directories(shared PUBLIC "${HDF5_INCLUDE_DIRS}")
target_include_directories(static PUBLIC "${HDF5_INCLUDE_DIRS}")

# Add binary directory to find generated version.h (cf cmake tutorial/Step1!) 
# TODO: Currently doesn't work but I don't know why ! :-(
#target_include_directories(shared PUBLIC "${PROJECT_BINARY_DIR}")
#target_include_directories(static PUBLIC "${PROJECT_BINARY_DIR}")

# Debug message
#get_target_property(debug shared INTERFACE_INCLUDE_DIRECTORIES)
#message("${debug}")

# Rename the output library names
set_target_properties(shared PROPERTIES OUTPUT_NAME gstlearn) # Will produce libgstlearn.so (on Linux)
set_target_properties(static PROPERTIES OUTPUT_NAME gstlearn) # Will produce libgstlearn.a  (on Linux)

# Set library version
set_target_properties(shared PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(static PROPERTIES VERSION ${PROJECT_VERSION})

# Anyone who wants to build a binary linking to gstlearn needs -fPIC
set_property(TARGET shared PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET static PROPERTY POSITION_INDEPENDENT_CODE 1)

# Link to boost
# Target for header-only dependencies. (Boost include directory)
# It should be PRIVATE if no headers of the gstlearn include boost files
target_link_libraries(shared PUBLIC Boost::boost)
target_link_libraries(static PUBLIC Boost::boost)

# Link to HDF5
target_link_libraries(shared PUBLIC ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(static PUBLIC ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES}) # Not a real link. Transmit info to those who link gstlearn.

############################## End loop on flavor

# Only compile shared library by default (exclude static from all)
set_target_properties(static PROPERTIES EXCLUDE_FROM_ALL TRUE)

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


# Install dynamic library in DESTINATION/lib folder
install(TARGETS shared EXPORT corelibs
        LIBRARY DESTINATION lib
        OPTIONAL
)
# Install static library in DESTINATION/lib folder
# TODO: After compiling static library only, 'make install' forces dynamic library reconstruction!
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

#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "gstlearn::shared" for configuration "Release"
set_property(TARGET gstlearn::shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(gstlearn::shared PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "NLopt::nlopt"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgstlearn.so.1.10.0a2"
  IMPORTED_SONAME_RELEASE "libgstlearn.so.1"
  )

list(APPEND _cmake_import_check_targets gstlearn::shared )
list(APPEND _cmake_import_check_files_for_gstlearn::shared "${_IMPORT_PREFIX}/lib/libgstlearn.so.1.10.0a2" )

# Import target "gstlearn::static" for configuration "Release"
set_property(TARGET gstlearn::static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(gstlearn::static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgstlearn.a"
  )

list(APPEND _cmake_import_check_targets gstlearn::static )
list(APPEND _cmake_import_check_files_for_gstlearn::static "${_IMPORT_PREFIX}/lib/libgstlearn.a" )

# Import target "gstlearn::gmtsph" for configuration "Release"
set_property(TARGET gstlearn::gmtsph APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(gstlearn::gmtsph PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgmtsph.a"
  )

list(APPEND _cmake_import_check_targets gstlearn::gmtsph )
list(APPEND _cmake_import_check_files_for_gstlearn::gmtsph "${_IMPORT_PREFIX}/lib/libgmtsph.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)

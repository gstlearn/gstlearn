# Default GNU installed directory names
include(GNUInstallDirs)

####################################################
## INSTALLATION

# Setup the installation directory
# https://stackoverflow.com/questions/39481958/setting-cmake-install-prefix-from-cmakelists-txt-file/39485990#39485990
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(DEFINED GSTLEARN_INSTALL_DIR) # Check for CMake variable first
    set(CMAKE_INSTALL_PREFIX ${GSTLEARN_INSTALL_DIR} CACHE PATH "" FORCE)
  else()
    if(DEFINED ENV{GSTLEARN_INSTALL_DIR}) # Check for environment variable then
      set(CMAKE_INSTALL_PREFIX $ENV{GSTLEARN_INSTALL_DIR} CACHE PATH "" FORCE)
    else()
      # Default installation folder
      set(CMAKE_INSTALL_PREFIX $ENV{HOME}/${PROJECT_NAME}_install CACHE PATH "" FORCE)
    endif()
  endif()
endif()

# Install the shared library in DESTINATION/lib folder
install(TARGETS shared
  EXPORT ${PROJECT_NAME}_corelibs
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
)

# Include directories
target_include_directories(shared PUBLIC
  # Installed includes are made PUBLIC for client who links the shared library
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

# Install the includes
install(
  DIRECTORY   ${INCLUDES}/            # Install library headers
  DESTINATION include/${PROJECT_NAME} # in DESTINATION/include/${PROJECT_NAME}
)

# Install the export file
install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_export.hpp
        DESTINATION include/${PROJECT_NAME}
)

# Install the version file
install(
  FILES ${PROJECT_BINARY_DIR}/version.h
  DESTINATION include/${PROJECT_NAME}
)

# Export the shared library cmake configuration (See above for corelibs definition)
install(
  EXPORT ${PROJECT_NAME}_corelibs
  FILE ${PROJECT_NAME}Config.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)

# Install doxygen html directories (optional)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen/html/                   # HTML files from build folder
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/doxygen/html # in DESTINATION/share
        OPTIONAL
)

####################################################
## UNINSTALLATION

# Add uninstall target
# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

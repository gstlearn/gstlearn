# Default GNU installed directory names
include(GNUInstallDirs)

include(CMakePackageConfigHelpers)

####################################################
## INSTALLATION

# Setup the installation directory
# https://stackoverflow.com/questions/39481958/setting-cmake-install-prefix-from-cmakelists-txt-file/39485990#39485990
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(DEFINED PROJECT_INSTALL_DIR) # Check for CMake variable first
    set(CMAKE_INSTALL_PREFIX ${PROJECT_INSTALL_DIR} CACHE PATH "" FORCE)
  else()
    if(DEFINED ENV{PROJECT_INSTALL_DIR}) # Check for environment variable then
      set(CMAKE_INSTALL_PREFIX $ENV{PROJECT_INSTALL_DIR} CACHE PATH "" FORCE)
    else()
      # Default installation folder
      set(CMAKE_INSTALL_PREFIX $ENV{HOME}/${PROJECT_NAME}_install CACHE PATH "" FORCE)
    endif()
  endif()
endif()

# Include directories
target_include_directories(shared PUBLIC
  # Installed includes are made PUBLIC for client who links the shared library
  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>"
  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}/${PROJECT_NAME}>"
)

# Install the shared library
install(
  TARGETS shared
  EXPORT ${PROJECT_NAME}_corelibs
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Install the includes
install(
  DIRECTORY   ${INCLUDES}
  DESTINATION include/${PROJECT_NAME}/${PROJECT_NAME}
)

# Install the export header file
install(FILES ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_export.hpp
        DESTINATION include/${PROJECT_NAME}/${PROJECT_NAME}
)

# Install the version header file
install(
  FILES ${PROJECT_BINARY_DIR}/version.h
  DESTINATION include/${PROJECT_NAME}/${PROJECT_NAME}
)

# Install doxygen html directories (optional)
install(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen/html/
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/doxygen/html
  OPTIONAL
)

# Export the shared library cmake configuration (See above for corelibs definition)
install(
  EXPORT ${PROJECT_NAME}_corelibs
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)

# Create the Config.cmake file
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)

# Install the Config.cmake file
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
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

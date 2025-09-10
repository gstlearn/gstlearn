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

target_include_directories(static PUBLIC
  # Installed includes are made PUBLIC for client who links the static library
  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>"
  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}/${PROJECT_NAME}>"
)

# Install the library.
#
# The static library is only installed if explicitely requested with:
#   cmake --install /path/to/build/dir --component static
# This is to make it possible to install it, but only if it has been compiled
# (to save time), though maybe that's trying too hard to be nice.
#
# This only works on Linux because on Windows the ARCHIVE part is used both
# for static libraries and dll export files (.lib for both). The latter must
# always be installed, so the former cannot be skipped, so the whole static
# target is skipped in this case (only the shared target can be installed).
# I do not know how to fix this.
#
# LIBRARY is for .so so not needed on Windows, RUNTIME is for .dll so not
# needed on not-Windows.
if (WIN32)
  install(
    TARGETS shared
    EXPORT ${PROJECT_NAME}_corelibs
    RUNTIME DESTINATION lib # for DLL, default destination is bin
    ARCHIVE DESTINATION lib # destination is the default (lib) but needed to install both
  )
else()
  # gmtsph is built by gstlearn itself (and always as a static library), so it
  # can be installed inside gstlearn itself. NLopt (if needed) is an external
  # dependency that users of gstlearn will have to provide themselves.
  install(
    TARGETS shared static gmtsph
    EXPORT ${PROJECT_NAME}_corelibs
    LIBRARY
    ARCHIVE COMPONENT static EXCLUDE_FROM_ALL
  )
endif()

# Install the includes
install(
  DIRECTORY   ${INCLUDES}
  DESTINATION include/${PROJECT_NAME}/${PROJECT_NAME}
  MESSAGE_NEVER
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

# Set whether the package will depend on NLopt.
set(NLOPT_REQUIRED TRUE)
get_target_property(nlopt_lib_type NLopt::nlopt TYPE)
if (${nlopt_lib_type} MATCHES STATIC_LIBRARY)
  # NLopt is still required for gstlearn::static but installing it is optional,
  # so ignore it. Users will get an error message for a missing dependency when
  # using gstlearn::static rather than a more helpful error message about the
  # NLopt package not being found when doing a find_package(gstlearn) but I
  # don't see how to make it better.
  set(NLOPT_REQUIRED FALSE)
endif()

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

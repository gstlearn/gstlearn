project(gstlearn CXX)

cmake_minimum_required(VERSION 3.10)

include(GNUInstallDirs)

# Add c++11 support whatever the compiler
# https://stackoverflow.com/questions/10851247/how-do-i-activate-c-11-in-cmake
set(CMAKE_CXX_STANDARD 11)

# Warning fiesta! (global property)
# https://cmake.org/cmake/help/latest/command/add_compile_options.html
if (MSVC)
  # Warning level 4 (4 = maximum, 0 = none)
  add_compile_options(/W4)
else()
  # Lots of warnings (-Wall = add some warnings, -Wextra = add a ton of warnings)
  add_compile_options(-Wall)
endif()

# Impose 'd' suffix in debug (global property)
set(CMAKE_DEBUG_POSTFIX d)

# Process sub directories
add_subdirectory(src)
add_subdirectory(doxygen)
add_subdirectory(tests)

# Export the package
install(
  EXPORT corelibs
  FILE gstlearnConfig.cmake
  NAMESPACE gstlearn::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/cmake
)
# Uncomment the following line if you want to register the package globally
# This will avoid the manual specification of gstlearn_DIR
# Register the gstlearn library
#set(CMAKE_EXPORT_PACKAGE_REGISTRY ON) # with CMake >= 3.15
#export(PACKAGE gstlearn)

# Uninstall target
# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()


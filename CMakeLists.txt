# Minimum CMake version
cmake_minimum_required(VERSION 3.19)

# Permits using project VERSION
cmake_policy(SET CMP0048 NEW)
project(gstlearn 
        VERSION      0.1.0 # Make it consistent with CHANGES file
        DESCRIPTION  "Geostatistics & Machine Learning toolbox"
        HOMEPAGE_URL https://github.com/gstlearn/gstlearn
        LANGUAGES    C CXX) # Enable C language for third party libraries
# Get the date in gstlearn_DATE variable
string(TIMESTAMP gstlearn_DATE "%Y/%m/%d")
# Fix version.h automatically
configure_file(version.h.in version.h)

# Default GNU installed directory names
include(GNUInstallDirs)

# Add test target (execution)
include(CTest)
enable_testing()

# Make Release version the default
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add c++11 support whatever the compiler
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Warning fiesta! (global property)
# https://cmake.org/cmake/help/latest/command/add_compile_options.html
if (MSVC)
  # Warning level 4 (4 = maximum, 0 = none)
  add_compile_options(/W4) 
else()
  # Lots of warnings (-Wall = add some warnings, -Wextra = add a ton of warnings)
  add_compile_options(-Wextra)
endif()

# Impose 'd' suffix in debug (global property)
set(CMAKE_DEBUG_POSTFIX d)

# Process sub directories
add_subdirectory(src)
add_subdirectory(doc)
add_subdirectory(tests)
add_subdirectory(doxygen)

# Export the package
# TODO : importing the package gstlearn imposes that static AND shared are installed
install(
  EXPORT corelibs # See src/CMakeLists.txt
  FILE gstlearnConfig.cmake
  NAMESPACE gstlearn::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/cmake
)

# Uncomment the following lines if you want to register the package globally
# This will avoid the manual specification of gstlearn_DIR
# Register the gstlearn library
# TODO: It registers the build directory - big deal... but latter find_package doesn't find gstlearn
#horlick:~> more .cmake/packages/gstlearn/0803e28e1916acff1b305242a59a11c4
#/home/fors/Projets/gstlearn/gstlearn/build
#cmake_policy(SET CMP0090 NEW)
#set(CMAKE_EXPORT_PACKAGE_REGISTRY ON) # with CMake >= 3.15
#export(PACKAGE gstlearn)

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

# Nothing is built by default (nothing in target 'all')
set(HELP_MESSAGE
"                                                             \r\n\
 Thank you for compiling gstlearn!                            \r\n\
 Please, use one of the available targets:                    \r\n\
  - shared     Build shared library                           \r\n\
  - static     Build static library                           \r\n\
  - doxygen    Build doxygen documentation                    \r\n\
  - build_test Build non-regression tests executables         \r\n\
  - test       Execute non-regression tests                   \r\n\
  - install    Install libraries, headers and documentation   \r\n\
  - uninstall  Uninstall libraries, headers and documentation \r\n\
  - clean      Clean generated files                          \r\n\
                                                              \r\n"
)

if (WIN32)
  message(${HELP_MESSAGE})
else()
  # To make output colored
  execute_process(COMMAND
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan ${HELP_MESSAGE}
    OUTPUT_FILE /dev/stdout ERROR_FILE /dev/stderr
  )
endif()

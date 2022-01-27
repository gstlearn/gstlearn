# For recent UseSWIG
cmake_minimum_required(VERSION 3.20)

# Permit using project VERSION
project(gstlearn
        # TODO : generate gstlearn.bib automatically
        VERSION      0.1.2 # Make it consistent with CHANGES and gstlearn.bib files
        DESCRIPTION  "Geostatistics & Machine Learning toolbox"
        HOMEPAGE_URL https://github.com/gstlearn/gstlearn
        LANGUAGES    C CXX) # Enable C language for third party libraries

####################################################
## CONFIGURATION

# Get the date in gstlearn_DATE variable
string(TIMESTAMP ${PROJECT_NAME}_DATE "%Y/%m/%d - %H:%M")

# Convert project name to uppercase
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UP)

# Fix version.h automatically
configure_file(version.h.in version.h)


####################################################
## BUILD & INSTALL

include(cmake/cpp.cmake)
include(cmake/doxygen.cmake)
include(cmake/doc.cmake)
include(cmake/install.cmake)


####################################################
## PACKAGES

add_subdirectory(python)


####################################################
## TESTING

# Add test target (execution)
include(CTest)
enable_testing()

add_subdirectory(tests)

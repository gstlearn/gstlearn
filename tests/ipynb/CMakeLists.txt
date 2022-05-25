if(NOT BUILD_PYTHON)
  return()
endif()

# TODO : All CMakeLists.txt for non-regression test scripts (py, R, ipynb) are similars (factorization ?)

# Look for Python3 (
find_package(Python3 REQUIRED)
get_filename_component(PYDIR ${Python3_EXECUTABLE} DIRECTORY)
set(PYTHON3_EXECUTABLE_DIR ${PYDIR} CACHE STRING "python3_dir")

set(CMAKE_FIND_DEBUG_MODE TRUE)
execute_process(COMMAND which jupyter-nbconvert)

# Look for jupyter-nbconvert executable
find_program(JUPYTER_NBCONVERT_EXECUTABLE REQUIRED
             NAMES jupyter-nbconvert
             PATHS ${PYTHON3_EXECUTABLE_DIR})
             #NO_DEFAULT_PATH)

# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil

# With GLOB:
file(GLOB TEST_SOURCES *.ipynb)
# Without GLOB:
#set(TEST_SOURCES 
#  DbNaming.ipynb)

# Generation folder
cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE}
           OUTPUT_VARIABLE GLOBAL_DST_DIR)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${GLOBAL_DST_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_DST_DIR})

# Create the output directory for logs
cmake_path(APPEND GLOBAL_DST_DIR output
           OUTPUT_VARIABLE TEST_DST_DIR)
file(MAKE_DIRECTORY ${TEST_DST_DIR})

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Run each jupyter notebook scripts
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Build the path to the jupyter notebook test script
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR ${TEST_SOURCE_FILE}
               OUTPUT_VARIABLE TEST_SOURCE_FILE)
    # Run the test (and generate asciidoc in ${TEST_DST_DIR} - asciidoc extension will be added by nbconvert) 
    add_test(NAME ${TEST_NAME}_ipynb
             COMMAND ${JUPYTER_NBCONVERT_EXECUTABLE} ${TEST_SOURCE_FILE} --to asciidoc --no-prompt --log-level=WARN --execute --output ${TEST_DST_DIR}/${TEST_NAME})
    # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
    # Use git diff (cross-platform) with no-index (no need to have files under git control) and ignore whitespaces (-w)
    add_test(NAME ${TEST_NAME}_ipynb_cmp
             COMMAND git diff --no-index -w ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}.asciidoc ${TEST_DST_DIR}/${TEST_NAME}.asciidoc)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# Look parent CMakeLists for MY_CTEST_COMMAND definition
add_custom_target(check_ipynb
                  COMMAND ${MY_CTEST_COMMAND})

# Add dependency for check_ipynb
add_dependencies(check_ipynb python_install)

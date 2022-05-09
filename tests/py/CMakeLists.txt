if(NOT BUILD_PYTHON)
  return()
endif()

# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil

# Look for Python3
find_package(Python3 REQUIRED)

# With GLOB:
file(GLOB TEST_SOURCES test*.py)
# Without GLOB:
#set(TEST_SOURCES 
#  test_Arguments.py)

# Generation folder
cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE}
           OUTPUT_VARIABLE GLOBAL_DST_DIR)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${GLOBAL_DST_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_DST_DIR})

# Create the output directory (logs)
cmake_path(APPEND GLOBAL_DST_DIR output
           OUTPUT_VARIABLE TEST_DST_DIR)
file(MAKE_DIRECTORY ${TEST_DST_DIR})

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Run each registered executable
set(TARGETS_SCRIPT "")
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Build the target test name
    set(TARGET_TEST_NAME ${TEST_NAME}_py)
    # Add to script targets list
    list(APPEND TARGETS_SCRIPT ${TARGET_TEST_NAME})
    # Run the test (and generate *.out in ${TEST_DST_DIR})
    add_test(NAME ${TARGET_TEST_NAME}
             COMMAND ${Python3_EXECUTABLE} ${TEST_SOURCE_FILE}
             WORKING_DIRECTORY ${TEST_DST_DIR})
    # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
    add_test(NAME ${TARGET_TEST_NAME}_cmp
             COMMAND git diff -w ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}.ref ${TEST_DST_DIR}/${TEST_NAME}.out
             WORKING_DIRECTORY ${TEST_DST_DIR})
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# Look parent CMakeLists for MY_CTEST_COMMAND definition
add_custom_target(check_py
  COMMAND ${MY_CTEST_COMMAND})

# Add dependency for check_py
add_dependencies(check_py python_install)

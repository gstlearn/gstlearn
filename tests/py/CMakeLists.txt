if(NOT BUILD_PYTHON)
  return()
endif()

# TODO : CMakeLists.txt for non-regression test scripts (py, R, ipynb, rmd) are similar (factorization ?)

# With GLOB:
# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
file(GLOB TEST_SOURCES test*.py)

# Without GLOB:
#set(TEST_SOURCES 
#    test_Arguments.py)

# Generation folder
if (IS_MULTI_CONFIG)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR $<CONFIG>
             OUTPUT_VARIABLE TEST_DST_DIR)
else()
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE}
             OUTPUT_VARIABLE TEST_DST_DIR)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_DST_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TEST_DST_DIR})

# Create the output directory for logs
add_custom_target(prepare_check_py
                  COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}")

# Add dependency for prepare_check_py
add_dependencies(prepare_check_py python_install)

# Create the main target
# Look parent CMakeLists for MY_CTEST_COMMAND definition
add_custom_target(check_py
                  COMMAND ${MY_CTEST_COMMAND})

# Add dependency for check_py
add_dependencies(check_py prepare_check_py)

# Look for Python3 (needed here because Python3_EXECUTABLE not populated everywhere)
find_package(Python3 REQUIRED)
# Python script test runner path
cmake_path(APPEND CMAKE_SOURCE_DIR python run_test_py.py
           OUTPUT_VARIABLE RUN_TEST_SCRIPT)

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Run each python scripts
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Build the path to the python test script
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR ${TEST_SOURCE_FILE}
               OUTPUT_VARIABLE TEST_SOURCE_FILE)
    # Run the test (and generate log in ${TEST_DST_DIR}) (use a test runner script for stdout redirection)
    add_test(NAME ${TEST_NAME}_py
             COMMAND ${Python3_EXECUTABLE} ${RUN_TEST_SCRIPT} ${TEST_SOURCE_FILE} ${TEST_DST_DIR})
    # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
    set(REF_FILE ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}.ref)
    # Use maybe a specific output log for the local computer
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR myoutput ${TEST_NAME}.ref
                OUTPUT_VARIABLE REF_FILE_LOCAL)
    if(EXISTS ${REF_FILE_LOCAL})
      set(REF_FILE ${REF_FILE_LOCAL})
    endif()
    # Use maybe a specific output log for MSVC
    if (MSVC)
      set(REF_FILE_MSVC ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}_msvc.ref)
      if(EXISTS ${REF_FILE_MSVC})
        set(REF_FILE ${REF_FILE_MSVC})
      endif()
    endif()
    # Use maybe a specific output log for clang
    if (CLANG)
      set(REF_FILE_CLANG ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}_clang.ref)
      if(EXISTS ${REF_FILE_CLANG})
        set(REF_FILE ${REF_FILE_CLANG})
      endif()
    endif()
    # Use git diff (cross-platform) with:
    # - no-index (no need to have files under git control)
    # - ignore lines starting with #NO_DIFF#
    # - ignore ALL whitespaces (https://stackoverflow.com/a/65462731)
    add_test(NAME ${TEST_NAME}_py_cmp
             COMMAND git diff --no-index -I "#NO_DIFF#" --ignore-blank-lines --ignore-all-space --ignore-cr-at-eol ${REF_FILE} ${TEST_DST_DIR}/${TEST_NAME}.out)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})


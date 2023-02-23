if(NOT BUILD_R)
  return()
endif()

# TODO : CMakeLists.txt for non-regression test scripts (py, R, ipynb, rmd) are similar (factorization ?)

# With GLOB:
# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
file(GLOB TEST_SOURCES test*.R)
# Without GLOB:
#set(TEST_SOURCES 
#    test_Arguments.R)

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

# Create the target and the output directory for logs
# Look parent CMakeLists for MY_CTEST_COMMAND definition
add_custom_target(check_r
                  COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
                  COMMAND ${MY_CTEST_COMMAND})

# Add dependency for check_r
add_dependencies(check_r r_install)

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Run each R scripts
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Build the path to the r test script
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR ${TEST_SOURCE_FILE}
               OUTPUT_VARIABLE TEST_SOURCE_FILE)
    # Run the test (and generate *.out in ${TEST_DST_DIR})
    add_test(NAME ${TEST_NAME}_r
             COMMAND ${R_EXECUTABLE} CMD BATCH --no-echo --no-timing ${TEST_SOURCE_FILE} ${TEST_DST_DIR}/${TEST_NAME}.out)
    # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
    # Use git diff (cross-platform) with no-index (no need to have files under git control) and ignore white spaces (-w)
    add_test(NAME ${TEST_NAME}_r_cmp
             COMMAND git diff --no-index -w ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}.ref ${TEST_DST_DIR}/${TEST_NAME}.out)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})







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

# Remove test_one
list(FILTER TEST_SOURCES EXCLUDE REGEX "test_one.R")

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
add_custom_target(prepare_check_r
                  COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}")

# Add dependency for prepare_check_r
add_dependencies(prepare_check_r r_install)

# Create the main target
# Look parent CMakeLists for MY_CTEST_COMMAND definition
add_custom_target(check_r
                  COMMAND ${MY_CTEST_COMMAND})

# Add dependency for check_r
add_dependencies(check_r prepare_check_r)

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
    set(REF_FILE ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}.ref)
    # Use maybe a specific output log for MinGW
    if (MINGW)
      set(REF_FILE_MSYS ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}_msys.ref)
      if(EXISTS ${REF_FILE_MSYS})
        set(REF_FILE ${REF_FILE_MSYS})
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
    add_test(NAME ${TEST_NAME}_r_cmp
             COMMAND git diff --no-index -I "#NO_DIFF#" --ignore-blank-lines --ignore-all-space --ignore-cr-at-eol ${REF_FILE} ${TEST_DST_DIR}/${TEST_NAME}.out)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})







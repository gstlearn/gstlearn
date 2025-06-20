if(NOT BUILD_R)
  return()
endif()

# TODO : CMakeLists.txt for non-regression test scripts (py, R, ipynb, rmd) are similar (factorization ?)
# R Markdown are located in doc/ folder
# They are not executed automatically for non-regression (see parent CMakeLists)

# With GLOB:
# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
cmake_path(APPEND PROJECT_SOURCE_DIR doc demo r
           OUTPUT_VARIABLE DEMO_DIR)
file(GLOB DEMO_SOURCES ${DEMO_DIR}/*.Rmd)
cmake_path(APPEND PROJECT_SOURCE_DIR doc courses r
           OUTPUT_VARIABLE COURSES_DIR)
file(GLOB COURSES_SOURCES ${COURSES_DIR}/*.Rmd)
set(TEST_SOURCES ${DEMO_SOURCES})
list(APPEND TEST_SOURCES ${COURSES_SOURCES})
# Without GLOB:
#set(TEST_SOURCES 
#    ../../doc/demo/r/Tuto_2D.Rmd)

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
add_custom_target(prepare_check_rmd
                  COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}")

# Create the target and the output directory for logs
# Look parent CMakeLists for MY_CTEST_COMMAND definition
add_custom_target(check_rmd
  COMMAND ${MY_CTEST_COMMAND}
  JOB_POOL console
)

# Add dependency for check_rmd
add_dependencies(check_rmd prepare_check_rmd r_install)

# Rmd script test runner path
cmake_path(APPEND CMAKE_SOURCE_DIR r run_test_rmd.R
           OUTPUT_VARIABLE RUN_TEST_SCRIPT)

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Run each R Markdown scripts
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Build the path to the R Markdown script
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR ${TEST_SOURCE_FILE}
               OUTPUT_VARIABLE TEST_SOURCE_FILE)
    # Run the test (and generate log in ${TEST_DST_DIR})
    add_test(NAME ${TEST_NAME}_rmd
             COMMAND ${R_EXECUTABLE} CMD BATCH --no-save --no-restore "--args ${TEST_SOURCE_FILE} ${TEST_DST_DIR}" ${RUN_TEST_SCRIPT} ${TEST_DST_DIR}/${TEST_NAME}.Rout)
    # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
    set(REF_FILE ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}.out)
    # Use maybe a specific output log for the local computer
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR myoutput ${TEST_NAME}.out
                OUTPUT_VARIABLE REF_FILE_LOCAL)
    if(EXISTS ${REF_FILE_LOCAL})
      set(REF_FILE ${REF_FILE_LOCAL})
    endif()
    # Use maybe a specific output log for MinGW
    if (MINGW)
      set(REF_FILE_MSYS ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}_msys.out)
      if(EXISTS ${REF_FILE_MSYS})
        set(REF_FILE ${REF_FILE_MSYS})
      endif()
    endif()
    # Use maybe a specific output log for clang
    if (CLANG)
      set(REF_FILE_CLANG ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}_clang.out)
      if(EXISTS ${REF_FILE_CLANG})
        set(REF_FILE ${REF_FILE_CLANG})
      endif()
    endif()
    # Use git diff (cross-platform) with:
    # - no-index (no need to have files under git control)
    # - ignore lines starting with #NO_DIFF#
    # - ignore ALL whitespaces (https://stackoverflow.com/a/65462731)
    add_test(NAME ${TEST_NAME}_rmd_cmp
             COMMAND git diff --no-index -I "#NO_DIFF#" --ignore-blank-lines --ignore-all-space --ignore-cr-at-eol ${REF_FILE} ${TEST_DST_DIR}/${TEST_NAME}.out)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

if(NOT BUILD_PYTHON)
  return()
endif()

# TODO : CMakeLists.txt for non-regression test scripts (py, R, ipynb, rmd) are similar (factorization ?)
# Jupyter notebooks are located in doc/ folder
# They are not executed automatically for non-regression (see parent CMakeLists)

# With GLOB:
# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
cmake_path(APPEND PROJECT_SOURCE_DIR doc demo python
           OUTPUT_VARIABLE DEMO_DIR)
file(GLOB DEMO_SOURCES ${DEMO_DIR}/*.ipynb)
cmake_path(APPEND PROJECT_SOURCE_DIR doc courses python
           OUTPUT_VARIABLE COURSES_DIR)
file(GLOB COURSES_SOURCES ${COURSES_DIR}/*.ipynb)
set(TEST_SOURCES ${DEMO_SOURCES})
list(APPEND TEST_SOURCES ${COURSES_SOURCES})
# Without GLOB:
#set(TEST_SOURCES
#    ../../doc/courses/python/05_Kriging.ipynb)

# Remove Tuto_SPDE and Tuto_SpatioTemp under MSVC (because scikits-sparse)
if (MSVC)
  list(FILTER TEST_SOURCES EXCLUDE REGEX ".*/Tuto_SPDE.ipynb")
  list(FILTER TEST_SOURCES EXCLUDE REGEX ".*/Tuto_SpatioTemp.ipynb")
endif()

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
add_custom_target(check_ipynb
                  COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
                  COMMAND ${MY_CTEST_COMMAND})

# Add dependency for check_ipynb
add_dependencies(check_ipynb python_install)

# Jupyter notebook script test runner path
cmake_path(APPEND CMAKE_SOURCE_DIR python run_test_ipynb.py
           OUTPUT_VARIABLE RUN_TEST_SCRIPT)

# Look for Python3 (needed here because Python3_EXECUTABLE not populated everywhere)
find_package(Python3 REQUIRED)
# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Run each jupyter notebook scripts
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Build the path to the jupyter notebook script
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR ${TEST_SOURCE_FILE}
               OUTPUT_VARIABLE TEST_SOURCE_FILE)
    # Run the test (and generate log in ${TEST_DST_DIR})
    add_test(NAME ${TEST_NAME}_ipynb
             COMMAND ${Python3_EXECUTABLE} ${RUN_TEST_SCRIPT} ${TEST_SOURCE_FILE} ${TEST_DST_DIR})
    # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
    set(REF_FILE ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}.asciidoc)
    # Use maybe a specific output log for MSVC
    if (MSVC)
      set(REF_FILE_MSVC ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}_msvc.asciidoc)
      if(EXISTS ${REF_FILE_MSVC})
        set(REF_FILE ${REF_FILE_MSVC})
      endif()
    endif()
    # Use maybe a specific output log for MSYS
    if (MINGW)
      set(REF_FILE_MSYS ${CMAKE_CURRENT_SOURCE_DIR}/output/${TEST_NAME}_msys.asciidoc)
      if(EXISTS ${REF_FILE_MSYS})
        set(REF_FILE ${REF_FILE_MSYS})
      endif()
    endif()
    # Use git diff (cross-platform) with:
    # - no-index (no need to have files under git control)
    # - ignore lines starting with #NO_DIFF#
    # - ignore whitespaces (-w) and blanck lines (-b) TODO : doesn't work ?
    add_test(NAME ${TEST_NAME}_ipynb_cmp
             COMMAND git diff --no-index -I "#NO_DIFF#" -b -w ${REF_FILE} ${TEST_DST_DIR}/${TEST_NAME}.asciidoc)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})


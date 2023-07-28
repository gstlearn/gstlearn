if(NOT BUILD_PYTHON)
  return()
endif()

# TODO : CMakeLists.txt for non-regression test scripts (py, R, ipynb, rmd) are similar (factorization ?)
# Jupyter notebooks are located in doc/ folder

# With GLOB:
# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
cmake_path(APPEND PROJECT_SOURCE_DIR doc courses python
           OUTPUT_VARIABLE COURSES_DIR)
file(GLOB COURSES_SOURCES ${COURSES_DIR}/*.ipynb)
set(COURSES_SOURCES ${COURSES_SOURCES})
# Without GLOB:
#set(COURSES_SOURCES
#    ../../doc/courses/python/01_gstlearn_start.ipynb)

# Generation folder
cmake_path(APPEND CMAKE_BINARY_DIR courses python
           OUTPUT_VARIABLE COURSES_DIR)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${COURSES_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${COURSES_DIR})

# Create the output directory for html
add_custom_target(prepare_build_courses_py
                  COMMAND ${CMAKE_COMMAND} -E make_directory "${COURSES_DIR}"
                  COMMENT "Creating ${COURSES_DIR}")
add_dependencies(prepare_build_courses_py python_install)

# Create build_courses_py target + add dependencies
add_custom_target(build_courses_py)
add_dependencies(build_courses_py prepare_build_courses_py)

# Jupyter notebook script runner path
cmake_path(APPEND CMAKE_SOURCE_DIR python run_test_ipynb.py
           OUTPUT_VARIABLE RUN_COURSES_SCRIPT)

# Look for Python3 (needed here because Python3_EXECUTABLE not populated everywhere)
find_package(Python3 REQUIRED)

# Run each jupyter notebook scripts
foreach(COURSES_SOURCE_FILE ${COURSES_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(COURSES_NAME ${COURSES_SOURCE_FILE} NAME_WE)
    # Build the path to the jupyter notebook script
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR ${COURSES_SOURCE_FILE}
               OUTPUT_VARIABLE COURSES_SOURCE_FILE)
    # Run the courses (and generate html in ${COURSES_DIR})
    add_custom_target(${COURSES_NAME}_ipynb_html
                      COMMAND ${Python3_EXECUTABLE} ${RUN_COURSES_SCRIPT} ${COURSES_SOURCE_FILE} ${COURSES_DIR} html)
                      
    # Add dependency
    add_dependencies(${COURSES_NAME}_ipynb_html prepare_build_courses_py)

    # Add dependency
    add_dependencies(build_courses_py ${COURSES_NAME}_ipynb_html)
    
endforeach(COURSES_SOURCE_FILE ${TEST_SOURCES})



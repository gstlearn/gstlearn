if(NOT BUILD_PYTHON)
  return()
endif()

# TODO : CMakeLists.txt for non-regression test scripts (py, R, ipynb, rmd) are similar (factorization ?)
# Jupyter notebooks are located in doc/ folder

# With GLOB:
# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
cmake_path(APPEND PROJECT_SOURCE_DIR doc demo python
           OUTPUT_VARIABLE DEMO_DIR)
file(GLOB DEMO_SOURCES ${DEMO_DIR}/*.ipynb)
set(DEMO_SOURCES ${DEMO_SOURCES})
# Without GLOB:
#set(DEMO_SOURCES
#    ../../doc/demo/python/Halieutic.ipynb)

# Generation folder
cmake_path(APPEND CMAKE_BINARY_DIR demo python
           OUTPUT_VARIABLE DEMO_DIR)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DEMO_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DEMO_DIR})

# Create the output directory for html
add_custom_target(prepare_build_demos_py
                  COMMAND ${CMAKE_COMMAND} -E make_directory "${DEMO_DIR}"
                  COMMENT "Creating ${DEMO_DIR}")
add_dependencies(prepare_build_demos_py python_install)

# Create build_demos_py target + add dependencies
add_custom_target(build_demos_py)
add_dependencies(build_demos_py prepare_build_demos_py)

# Jupyter notebook script runner path
cmake_path(APPEND CMAKE_SOURCE_DIR python run_test_ipynb.py
           OUTPUT_VARIABLE RUN_DEMO_SCRIPT)

# Look for Python3 (needed here because Python3_EXECUTABLE not populated everywhere)
find_package(Python3 REQUIRED)

# Run each jupyter notebook scripts
foreach(DEMO_SOURCE_FILE ${DEMO_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(DEMO_NAME ${DEMO_SOURCE_FILE} NAME_WE)
    # Build the path to the jupyter notebook script
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR ${DEMO_SOURCE_FILE}
               OUTPUT_VARIABLE DEMO_SOURCE_FILE)
    # Run the demo (and generate html in ${DEMO_DIR})
    add_custom_target(${DEMO_NAME}_ipynb_html
                      COMMAND ${Python3_EXECUTABLE} -Xfrozen_modules=off ${RUN_DEMO_SCRIPT} ${DEMO_SOURCE_FILE} ${DEMO_DIR} html)

    # Add dependency
    add_dependencies(${DEMO_NAME}_ipynb_html prepare_build_demos_py)

    # Add dependency
    add_dependencies(build_demos_py ${DEMO_NAME}_ipynb_html)
    
endforeach(DEMO_SOURCE_FILE ${TEST_SOURCES})



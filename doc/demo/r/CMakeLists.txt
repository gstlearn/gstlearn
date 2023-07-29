if(NOT BUILD_R)
  return()
endif()

# TODO : CMakeLists.txt for non-regression test scripts (py, R, ipynb, rmd) are similar (factorization ?)
# R Markdown are located in doc/ folder

# With GLOB:
# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
cmake_path(APPEND PROJECT_SOURCE_DIR doc demo r
           OUTPUT_VARIABLE DEMO_DIR)
file(GLOB DEMO_SOURCES ${DEMO_DIR}/*.Rmd)
set(DEMO_SOURCES ${DEMO_SOURCES})
# Without GLOB:
#set(DEMO_SOURCES
#    ../../doc/demo/r/Tuto_2D.Rmd)

# Generation folder
cmake_path(APPEND CMAKE_BINARY_DIR demo r
           OUTPUT_VARIABLE DEMO_DIR)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DEMO_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DEMO_DIR})

# Create the output directory for html
add_custom_target(prepare_build_demos_r
                  COMMAND ${CMAKE_COMMAND} -E make_directory "${DEMO_DIR}"
                  COMMENT "Creating ${DEMO_DIR}")
add_dependencies(prepare_build_demos_r r_install)

# Create build_demos_r target + add dependencies
add_custom_target(build_demos_r)

# R Markdown script runner path
cmake_path(APPEND CMAKE_SOURCE_DIR r run_test_rmd.R
           OUTPUT_VARIABLE RUN_DEMO_SCRIPT)

# Run each R Markdown scripts
foreach(DEMO_SOURCE_FILE ${DEMO_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(DEMO_NAME ${DEMO_SOURCE_FILE} NAME_WE)
    # Build the path to the R Markdown script
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR ${DEMO_SOURCE_FILE}
               OUTPUT_VARIABLE DEMO_SOURCE_FILE)
    # Run the demo (and generate html in ${DEMO_DIR})
    # https://stackoverflow.com/a/30554504/3952924
    add_custom_target(${DEMO_NAME}_rmd_html
                      COMMAND mkfifo ${DEMO_NAME}_fifo
                      COMMAND cat ${DEMO_NAME}_fifo &
                      COMMAND ${R_EXECUTABLE} CMD BATCH --no-save --no-restore "--args ${DEMO_SOURCE_FILE} ${DEMO_DIR} html" ${RUN_DEMO_SCRIPT} ${DEMO_NAME}_fifo
                      COMMAND rm ${DEMO_NAME}_fifo)

    # Add dependency for the Rmd
    add_dependencies(${DEMO_NAME}_rmd_html prepare_build_demos_r)
    
    # Add the Rmd to global dependency
    add_dependencies(build_demos_r ${DEMO_NAME}_rmd_html)
    
endforeach(DEMO_SOURCE_FILE ${TEST_SOURCES})



if(NOT BUILD_R)
  return()
endif()

# TODO : CMakeLists.txt for non-regression test scripts (py, R, ipynb, rmd) are similar (factorization ?)
# R Markdown are located in doc/ folder

# With GLOB:
# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
cmake_path(APPEND PROJECT_SOURCE_DIR doc courses r
           OUTPUT_VARIABLE COURSES_DIR)
file(GLOB COURSES_SOURCES ${COURSES_DIR}/*.Rmd)
set(COURSES_SOURCES ${COURSES_SOURCES})
# Without GLOB:
#set(COURSES_SOURCES
#    ../../doc/courses/r/00_Learn_R.Rmd)

# Generation folder
cmake_path(APPEND CMAKE_BINARY_DIR courses r
           OUTPUT_VARIABLE COURSES_DIR)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${COURSES_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${COURSES_DIR})

# Create the output directory for html
add_custom_target(prepare_build_courses_r
                  COMMAND ${CMAKE_COMMAND} -E make_directory "${COURSES_DIR}"
                  COMMENT "Creating ${COURSES_DIR}")
add_dependencies(prepare_build_courses_r r_install)

# Create bcourse_courses_r target + add dependencies
add_custom_target(build_courses_r)

# R Markdown script runner path
cmake_path(APPEND CMAKE_SOURCE_DIR r run_test_rmd.R
           OUTPUT_VARIABLE RUN_COURSES_SCRIPT)

# Run each R Markdown scripts
foreach(COURSES_SOURCE_FILE ${COURSES_SOURCES})
    # Retrieve script file name without extension
    get_filename_component(COURSES_NAME ${COURSES_SOURCE_FILE} NAME_WE)
    # Build the path to the R Markdown script
    cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR ${COURSES_SOURCE_FILE}
               OUTPUT_VARIABLE COURSES_SOURCE_FILE)
    # Run the courses (and generate html in ${COURSES_DIR})
    # https://stackoverflow.com/a/30554504/3952924
    add_custom_target(${COURSES_NAME}_rmd_html
                      COMMAND mkfifo ${COURSES_NAME}_fifo
                      COMMAND cat ${COURSES_NAME}_fifo &
                      COMMAND ${R_EXECUTABLE} CMD BATCH --no-save --no-restore "--args ${COURSES_SOURCE_FILE} ${COURSES_DIR} html" ${RUN_COURSES_SCRIPT} ${COURSES_NAME}_fifo
                      COMMAND rm ${COURSES_NAME}_fifo)

    # Add dependency for the Rmd
    add_dependencies(${COURSES_NAME}_rmd_html prepare_build_courses_r)
    
    # Add the Rmd to global dependency
    add_dependencies(build_courses_r ${COURSES_NAME}_rmd_html)
    
endforeach(COURSES_SOURCE_FILE ${TEST_SOURCES})



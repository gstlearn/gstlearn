if(CMAKE_CONFIGURATION_TYPES)
  # Currently these tests cannot run under Visual
  # TODO : make data non-reg tests available under Visual ($<CONFIG> for output directoryu)
  return()
endif()

# List special non-regression test sources
set(TEST_SOURCES 
  testim.cpp
  testpgs.cpp
  testvar.cpp)

# Generation folder (all tests in one folder)
if (CMAKE_CONFIGURATION_TYPES)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR $<CONFIG>
             OUTPUT_VARIABLE GLOBAL_DST_DIR)
else()
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE}
             OUTPUT_VARIABLE GLOBAL_DST_DIR)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${GLOBAL_DST_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_DST_DIR})

# Compile each test executable
set(TARGETS_EXE "")
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve source file name without extension (will become executable name)
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Add to executable targets list
    list(APPEND TARGETS_EXE ${TEST_NAME})
    # Define sources list for the target executable
    add_executable(${TEST_NAME} EXCLUDE_FROM_ALL ${TEST_SOURCE_FILE})
    # Link each test to shared library
    target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAME}::shared)
    # Trigger the build of the test with the target build_tests
    add_dependencies(build_tests ${TEST_NAME})
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Iterate over each target executables
foreach(TARGET_EXE ${TARGETS_EXE})
    set(BDIR "")
    if(${TARGET_EXE} MATCHES "testim*")
        set(BDIR "Jeu")
    endif()
    if(${TARGET_EXE} MATCHES "testpgs*")
        set(BDIR "PGS")
    endif()
    if(${TARGET_EXE} MATCHES "testvar*")
        set(BDIR "Var")
    endif()
    # Calculate full path of target executable
    cmake_path(APPEND GLOBAL_DST_DIR ${TARGET_EXE}
               OUTPUT_VARIABLE TARGET_EXE)
    # Retrieve all test cases
    file(GLOB TEST_SRC_DIRS ${BDIR}*)
    foreach(TEST_SRC_DIR ${TEST_SRC_DIRS})
        # Drop absolute directory path
        get_filename_component(DIR_NAME "${TEST_SRC_DIR}" NAME)
        # Calculate output directory
        cmake_path(APPEND GLOBAL_DST_DIR ${DIR_NAME}
                   OUTPUT_VARIABLE TEST_DST_DIR)
        file(MAKE_DIRECTORY ${TEST_DST_DIR}) # This contains a generator expression in multiconfig => Not possible
        # Run the test (and generate Result.out in ${TEST_DST_DIR})
        add_test(NAME ${DIR_NAME}
                 COMMAND ${TARGET_EXE} ${TEST_SRC_DIR}
                 WORKING_DIRECTORY ${TEST_DST_DIR})
        # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
        add_test(NAME ${DIR_NAME}_cmp
                 COMMAND git diff -w ${TEST_SRC_DIR}/Result.ref ${TEST_DST_DIR}/Result.out
                 WORKING_DIRECTORY ${TEST_DST_DIR})
    endforeach(TEST_SRC_DIR ${TEST_SRC_DIRS})
endforeach(TARGET_EXE ${TARGETS_EXE})

 # Look parent CMakeLists for MY_CTEST_COMMAND definition
if (WIN32)
  # Need to copy C++ shared library to tests directory
  add_custom_target(check_data
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}::shared> ${GLOBAL_DST_DIR}
    COMMAND ${MY_CTEST_COMMAND} DEPENDS ${TARGETS_EXE}
  )
else()
  add_custom_target(check_data
    COMMAND ${MY_CTEST_COMMAND} DEPENDS ${TARGETS_EXE}
  )
endif()

# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil

# With GLOB:
file(GLOB TEST_SOURCES test_*.cpp)
# Without GLOB:
#set(TEST_SOURCES 
#  test_Francky.cpp)
#  test_SPDEDrift.cpp)

# Add test_HDF5 if supported
if (USE_HDF5)
  list(APPEND TEST_SOURCES ztest_HDF5.cpp)
endif()

# Generation folder
if (CMAKE_CONFIGURATION_TYPES)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR $<CONFIG>
             OUTPUT_VARIABLE TEST_DST_DIR)
else()
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE}
             OUTPUT_VARIABLE TEST_DST_DIR)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_DST_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TEST_DST_DIR})

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

# Create the target and the output directory for logs
# Look parent CMakeLists for MY_CTEST_COMMAND definition
if (WIN32)
  # Need to copy C++ shared library to tests directory
  add_custom_target(check_cpp
    COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}::shared> ${TEST_DST_DIR}
    COMMAND ${MY_CTEST_COMMAND} DEPENDS ${TARGETS_EXE}
  )
else()
  add_custom_target(check_cpp
    COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
    COMMAND ${MY_CTEST_COMMAND} DEPENDS ${TARGETS_EXE}
  )
endif()

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Deactivate a particular test (which is under construction)
#list(FILTER TARGETS_EXE EXCLUDE REGEX "test_kd.*")

# Run each registered executable
foreach(TARGET_EXE ${TARGETS_EXE})
    # Run the test (and generate *.out in ${TEST_DST_DIR})
    add_test(NAME ${TARGET_EXE}
             COMMAND ${TARGET_EXE}
             WORKING_DIRECTORY ${TEST_DST_DIR})
    # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
    # Use git diff (cross-platform) with no-index (no need to have files under git control) and ignore whitespaces (-w)
    add_test(NAME ${TARGET_EXE}_cmp
             COMMAND git diff --no-index -w ${CMAKE_CURRENT_SOURCE_DIR}/output/${TARGET_EXE}.ref ${TEST_DST_DIR}/${TARGET_EXE}.out)
endforeach(TARGET_EXE ${TARGETS_EXE})



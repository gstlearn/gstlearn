# List the source files
set(TEST_SOURCES 
  testInter_Ask.cpp
  testInter_GibbsMMulti.cpp)

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

# Compile each test executable
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve source file name without extension (will become executable name)
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Define sources list for the target executable
    add_executable(${TEST_NAME} ${TEST_SOURCE_FILE})
    # Link gstlearn to each test
    target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAME}::shared)
    # Trigger the build of the test with the target build_tests
    add_dependencies(build_tests ${TEST_NAME})
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# No automatic tests for interactive executable

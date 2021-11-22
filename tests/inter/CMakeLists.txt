# List the source files
set(TEST_SOURCES 
  testInter_Ask.cpp
  testInter_GibbsMMulti.cpp)

# Compile each executable
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve source file name without extension (will become executable name)
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Define sources list for the target executable
    add_executable(${TEST_NAME} ${TEST_SOURCE_FILE})
    # Include directories (-I)
    target_include_directories(${TEST_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>")
    # Make sure gstlearn is linked to each test
    target_link_libraries(${TEST_NAME} shared)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# No automatic tests for interactive executable
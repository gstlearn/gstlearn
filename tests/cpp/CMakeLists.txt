# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil

# With GLOB:
file(GLOB TEST_SOURCES test*.cpp)
# Without GLOB:
#set(TEST_SOURCES test_AllGibbs.cpp test_krige.cpp) # To be completed...

# Compile each executable
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve source file name without extension (will become executable name)
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Append .exe extension
    set(TEST_NAME ${TEST_NAME}.exe)
    # Define sources list for the target executable
    add_executable(${TEST_NAME} ${TEST_SOURCE_FILE})
    # Include directories (-I)
    target_include_directories(${TEST_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>")
    # Make sure gstlearn is linked to each test
    target_link_libraries(${TEST_NAME} gstlearn-shared)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# Display test output in case of failure
set (CTEST_OUTPUT_ON_FAILURE ON)

# With GLOB:
file(GLOB TEST_SOURCES test_*.cpp)  # Run only tests starting by 'test_'
# Without GLOB:
#set(TEST_SOURCES test_AllGibbs.cpp test_krige.cpp) # To be completed...

# Run each executable (starting by 'test_'
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve source file name without extension (will become executable name)
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
    # Append .exe extension
    set(TEST_NAME ${TEST_NAME}.exe)
    # Add current test (to be run)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})


# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil

# With GLOB:
file(GLOB TEST_SOURCES test*.cpp)
# Without GLOB:
#set(TEST_SOURCES test_AllGibbs.cpp test_krige.cpp) # To be completed...

# Prepare executable target list
set(TARGETS_EXE "")

# Compile each executable
foreach(TEST_SOURCE_FILE ${TEST_SOURCES})
    # Retrieve source file name without extension (will become executable name)
    get_filename_component(TEST_NAME ${TEST_SOURCE_FILE} NAME_WE)
	# Add to executable targets list
	list(APPEND TARGETS_EXE ${TEST_NAME})
    # Define sources list for the target executable
    add_executable(${TEST_NAME} ${TEST_SOURCE_FILE})
    # Include directories (-I)
    target_include_directories(${TEST_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>")
    # Make sure gstlearn is linked to each test
    target_link_libraries(${TEST_NAME} gstlearn-shared)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# Display test output in case of failure
set (CTEST_OUTPUT_ON_FAILURE ON)

# Run each registered executable
foreach(TARGET_EXE ${TARGETS_EXE})
    # Add current test (to be run)
    add_test(NAME ${TARGET_EXE} COMMAND ${TARGET_EXE})
endforeach(TARGET_EXE ${TARGETS_EXE})

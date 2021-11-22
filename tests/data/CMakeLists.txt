# List special non-regression test sources
set(TEST_SOURCES 
  testim.cpp
  testpgs.cpp
  testvar.cpp)

# Prepare executable targets list
set(TARGETS_EXE "")

# For each test executable
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
    target_link_libraries(${TEST_NAME} shared)
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

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
    # Retrieve all test cases
    file(GLOB TEST_DIRS ${BDIR}*)
    foreach(TEST_DIR ${TEST_DIRS})
        # Drop absolute directory path
        get_filename_component(DIR_NAME "${TEST_DIR}" NAME)
        # Run the test (The executable will generate Result.out in ${DIR_NAME})
        add_test(NAME ${DIR_NAME}
                 COMMAND ${TARGET_EXE} ${TEST_DIR}
                 WORKING_DIRECTORY ${TEST_DIR})
        # Compare the output result
        add_test(NAME ${DIR_NAME}_cmp
                 COMMAND ${CMAKE_COMMAND} -E compare_files ${TEST_DIR}/Result.ref ${TEST_DIR}/Result.out
                 WORKING_DIRECTORY ${TEST_DIR})
    endforeach(TEST_DIR ${TEST_DIRS})
endforeach(TARGET_EXE ${TARGETS_EXE})



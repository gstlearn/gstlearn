# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil

# With GLOB:
file(GLOB TEST_SOURCES test_*.cpp)
# Without GLOB:
#set(TEST_SOURCES 
#  test_Anam.cpp)
#  test_SPDEDrift.cpp)

# Remove test_ppmt (not yet functional)
#list(FILTER TEST_SOURCES EXCLUDE REGEX "test_ppmt.cpp")

# Remove test_HDF5 if not supported
if (NOT USE_HDF5)
  list(FILTER TEST_SOURCES EXCLUDE REGEX "test_HDF5.cpp")
endif()

# Exclude one test from output comparison
set(EXCLUDE_FROM_DIFF "test_a_template")

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

# Define the prepare target to create the output directory for logs
if (WIN32)
  # Need to copy C++ shared library to tests directory
  add_custom_target(prepare_cpp
    COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}::shared> ${TEST_DST_DIR}
  )
else()
  add_custom_target(prepare_cpp
    COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
  )
endif()

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
    # Trigger the prepare target each time a test is compiled
    add_dependencies(${TEST_NAME} prepare_cpp)
    # Trigger the build of the test with the target build_tests
    add_dependencies(build_tests ${TEST_NAME})
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Deactivate a particular test (which is under construction)
list(FILTER TARGETS_EXE EXCLUDE REGEX "test_a_template")

# Run each registered executable
foreach(TARGET_EXE ${TARGETS_EXE})
    # Run the test (and generate *.out in ${TEST_DST_DIR})
    add_test(NAME ${TARGET_EXE}
             COMMAND ${TARGET_EXE}
             WORKING_DIRECTORY ${TEST_DST_DIR})
    if (NOT ${TARGET_EXE} STREQUAL ${EXCLUDE_FROM_DIFF})
      # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
      set(REF_FILE ${CMAKE_CURRENT_SOURCE_DIR}/output/${TARGET_EXE}.ref)
      # Use maybe a specific output log for MSVC
      if (MSVC)
        set(REF_FILE_MSVC ${CMAKE_CURRENT_SOURCE_DIR}/output/${TARGET_EXE}_msvc.ref)
        if(EXISTS ${REF_FILE_MSVC})
          set(REF_FILE ${REF_FILE_MSVC})
        endif()
      endif()
      # Use maybe a specific output log for MinGW
      if (MINGW)
        set(REF_FILE_MSYS ${CMAKE_CURRENT_SOURCE_DIR}/output/${TARGET_EXE}_msys.ref)
        if(EXISTS ${REF_FILE_MSYS})
          set(REF_FILE ${REF_FILE_MSYS})
        endif()
      endif()
      # Use git diff (cross-platform) with no-index (no need to have files under git control) and ignore whitespaces (-w)
      add_test(NAME ${TARGET_EXE}_cmp
               COMMAND git diff --no-index -w  ${REF_FILE} ${TEST_DST_DIR}/${TARGET_EXE}.out)
    endif()
endforeach(TARGET_EXE ${TARGETS_EXE})

# Create the check target to launch the tests
# Look parent CMakeLists for MY_CTEST_COMMAND definition
add_custom_target(check_cpp
  COMMAND ${MY_CTEST_COMMAND} DEPENDS ${TARGETS_EXE}
)

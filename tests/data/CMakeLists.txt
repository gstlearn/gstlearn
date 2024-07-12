# List special non-regression test sources
set(TEST_SOURCES 
  testim.cpp
  testpgs.cpp
  testvar.cpp)

# Generation folder (all tests in one folder)
if (IS_MULTI_CONFIG)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR $<CONFIG>
             OUTPUT_VARIABLE GLOBAL_DST_DIR)
else()
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${CMAKE_BUILD_TYPE}
             OUTPUT_VARIABLE GLOBAL_DST_DIR)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${GLOBAL_DST_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_DST_DIR})

# Define the prepare target to create the output directory for logs
if (WIN32)
  # Need to copy C++ shared library to tests directory
  if(MSYS)
    add_custom_target(prepare_data
      COMMAND ${CMAKE_COMMAND} -E make_directory "${GLOBAL_DST_DIR}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}::shared> ${GLOBAL_DST_DIR}
      # have both lib${PROJECT_NAME}.dll and ${PROJECT_NAME}.dll in ${GLOBAL_DST_DIR}
      # this fixes a bug on MSYS when building R wrappers with shared
      # objects before linking C++ test binaries
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${GLOBAL_DST_DIR}/lib${PROJECT_NAME}.dll ${GLOBAL_DST_DIR}/${PROJECT_NAME}.dll
    )
  else()
    add_custom_target(prepare_data
      COMMAND ${CMAKE_COMMAND} -E make_directory "${GLOBAL_DST_DIR}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${PROJECT_NAME}::shared> ${GLOBAL_DST_DIR}
    )
  endif()
else()
  add_custom_target(prepare_data
    COMMAND ${CMAKE_COMMAND} -E make_directory "${GLOBAL_DST_DIR}"
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
    add_executable(${TEST_NAME} ${TEST_SOURCE_FILE})
    # Link each test to shared library
    target_link_libraries(${TEST_NAME} PRIVATE ${PROJECT_NAME}::shared)
    # Trigger the prepare target each time a test is compiled
    add_dependencies(${TEST_NAME} prepare_data)
    # Trigger the build of the test with the target build_tests
    add_dependencies(build_tests ${TEST_NAME})
endforeach(TEST_SOURCE_FILE ${TEST_SOURCES})

# Display test output in case of failure
set(CTEST_OUTPUT_ON_FAILURE ON)

# Iterate over each target executables
set(TARGETS_PREPARE "")
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
        # Define the prepare target to create the output directory for logs
        add_custom_target(${DIR_NAME}_prepare
          COMMAND ${CMAKE_COMMAND} -E make_directory "${TEST_DST_DIR}"
        )
        # Add to prepare targets list
        list(APPEND TARGETS_PREPARE ${DIR_NAME}_prepare)
        # Run the test (and generate Result.out in ${TEST_DST_DIR})
        add_test(NAME ${DIR_NAME}
                 COMMAND ${TARGET_EXE} ${TEST_SRC_DIR}
                 WORKING_DIRECTORY ${TEST_DST_DIR})
        # Compare the output result (diff output is shown when CTEST_OUTPUT_ON_FAILURE=1)
        set(REF_FILE ${TEST_SRC_DIR}/Result.ref)
        # Use maybe a specific output log for the local computer
        cmake_path(APPEND CMAKE_CURRENT_SOURCE_DIR myoutput ${DIR_NAME} Result.ref
                   OUTPUT_VARIABLE REF_FILE_LOCAL)
        if(EXISTS ${REF_FILE_LOCAL})
          set(REF_FILE ${REF_FILE_LOCAL})
        endif()
        # Use maybe a specific output log for MSVC
        if (MSVC)
          set(REF_FILE_MSVC ${TEST_SRC_DIR}/Result_msvc.ref)
          if(EXISTS ${REF_FILE_MSVC})
            set(REF_FILE ${REF_FILE_MSVC})
          endif()
        endif()
        # Use maybe a specific output log for MinGW
        if (MINGW)
          set(REF_FILE_MSYS ${TEST_SRC_DIR}/Result_msys.ref)
          if(EXISTS ${REF_FILE_MSYS})
            set(REF_FILE ${REF_FILE_MSYS})
          endif()
        endif()
        # Use maybe a specific output log for clang
        if (CLANG)
          set(REF_FILE_CLANG ${TEST_SRC_DIR}/Result_clang.ref)
          if(EXISTS ${REF_FILE_CLANG})
            set(REF_FILE ${REF_FILE_CLANG})
          endif()
        endif()
        # Use git diff (cross-platform) with:
        # - no-index (no need to have files under git control)
        # - ignore lines starting with #NO_DIFF#
        # - ignore ALL whitespaces (https://stackoverflow.com/a/65462731)
        add_test(NAME ${DIR_NAME}_cmp
                 COMMAND git diff --no-index -I "#NO_DIFF#" --ignore-blank-lines --ignore-all-space --ignore-cr-at-eol ${REF_FILE} ${TEST_DST_DIR}/Result.out
                 WORKING_DIRECTORY ${TEST_DST_DIR})
    endforeach(TEST_SRC_DIR ${TEST_SRC_DIRS})
endforeach(TARGET_EXE ${TARGETS_EXE})

 # Look parent CMakeLists for MY_CTEST_COMMAND definition
add_custom_target(check_data
  COMMAND ${MY_CTEST_COMMAND} DEPENDS ${TARGETS_PREPARE} ${TARGETS_EXE}
)

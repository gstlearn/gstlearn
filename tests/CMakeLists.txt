# Add 'build_tests' target for building all tests executable
add_custom_target(build_tests)

# Add 'check' target for building then running the tests
add_custom_target(check)

# Define ctest command # TODO : really needed ? (yaml use --target check)
if (MSVC) # TODO use IS_MULTI_CONFIG
  # Need configuration to be tested
#  set(MY_CTEST_COMMAND ${CMAKE_CTEST_COMMAND} --verbose -C $<CONFIG> --output-on-failure)
  set(MY_CTEST_COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIG> --output-on-failure)
else()
#  set(MY_CTEST_COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
  set(MY_CTEST_COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure)
endif()

# TODO : Output logs of each target are mixed when using several processors

add_subdirectory(data)
add_subdirectory(cpp)
add_subdirectory(inter)
set(CHECK_DEP check_cpp check_data)

if (BUILD_PYTHON)
  add_subdirectory(py)
  #add_subdirectory(ipynb)
  list(APPEND CHECK_DEP check_py)
endif()

if (BUILD_R)
  add_subdirectory(r)
  list(APPEND CHECK_DEP check_r)
endif()

add_dependencies(check ${CHECK_DEP})

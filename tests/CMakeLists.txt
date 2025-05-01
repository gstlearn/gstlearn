# Add 'build_tests' target for building all tests executable
add_custom_target(build_tests)

# Add 'check' target for building then running the tests
add_custom_target(check)

# Define ctest command
#set(VERBOSE_TEST --verbose)
if (IS_MULTI_CONFIG)
  set(CONFIG_TEST -C $<CONFIG>)
endif()
set(MY_CTEST_COMMAND ${CMAKE_CTEST_COMMAND} ${VERBOSE_TEST} ${CONFIG_TEST} --output-on-failure)

add_subdirectory(data)
add_subdirectory(cpp)
add_subdirectory(inter)
set(CHECK_DEP check_cpp check_data)

if (BUILD_PYTHON)
  add_subdirectory(py)
  list(APPEND CHECK_DEP check_py)
  add_subdirectory(ipynb)
  #list(APPEND CHECK_DEP check_ipynb)
endif()

if (BUILD_R)
  add_subdirectory(r)
  list(APPEND CHECK_DEP check_r)
  add_subdirectory(rmd)
  #list(APPEND CHECK_DEP check_rmd)
endif()

add_dependencies(check ${CHECK_DEP})

if(NOT BUILD_TESTS)
  return()
endif()

# Target for building the tests executable
add_custom_target(build_test ALL)

add_subdirectory(data)
add_subdirectory(cpp)
add_subdirectory(inter)

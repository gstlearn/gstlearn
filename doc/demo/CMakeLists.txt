# Add 'build_demos' target for processing all demo scripts
add_custom_target(build_demos)

set(BUILD_DEP "")

if (BUILD_PYTHON)
  add_subdirectory(python)
  list(APPEND BUILD_DEP build_demos_py)
endif()

if (BUILD_R)
  add_subdirectory(r)
  list(APPEND BUILD_DEP build_demos_r)
endif()

add_dependencies(build_demos ${BUILD_DEP})



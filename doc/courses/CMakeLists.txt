# Add 'build_courses' target for processing all courses scripts
add_custom_target(build_courses)

set(BUILD_DEP "")

if (BUILD_PYTHON)
  add_subdirectory(python)
  list(APPEND BUILD_DEP build_courses_py)
endif()

if (BUILD_R)
  add_subdirectory(r)
  list(APPEND BUILD_DEP build_courses_r)
endif()

if (BUILD_PYTHON OR BUILD_R)
  add_dependencies(build_courses ${BUILD_DEP})
endif()



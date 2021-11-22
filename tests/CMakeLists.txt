# Build the tests executable except if the user excludes them
if (NOT EXCLUDE_TEST)
  add_subdirectory(data)
  add_subdirectory(cpp)
  add_subdirectory(inter)
endif()

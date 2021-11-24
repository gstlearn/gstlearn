# Copy the data used by tests
if (NOT EXCLUDE_TEST)
  file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

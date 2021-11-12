# Using GLOB could be a bad idea (but use it for the tests)
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
file(GLOB TEST_SOURCES cpp/test*.cpp)
foreach(testsourcefile ${TEST_SOURCES})
    # Retrieve source file name without extension (will become executable name)
    get_filename_component(testname ${testsourcefile} NAME_WE)
    # Define sources list for the target executable
    add_executable(${testname} ${testsourcefile})
    # Include directories (-I)
    target_include_directories(${testname} PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>")
    # Make sure gstlearn is linked to each test
    target_link_libraries(${testname} gstlearn-shared)
endforeach(testsourcefile ${TEST_SOURCES})

# Without GLOB:
#set(TESTS_SRC Inter_Ask _AllGibbs _krige) # A completer...
#foreach(target ${TESTS_SRC})
#    add_executable(test${target} cpp/test${target}.cpp)
#    target_link_libraries(test${target} gstlearn-shared)
#endforeach()

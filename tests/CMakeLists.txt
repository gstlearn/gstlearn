# Using GLOB could be a bad idea

file(GLOB TEST_SOURCES cpp/test*.cpp)
foreach(testsourcefile ${TEST_SOURCES})
    #message("Processing " ${testsourcefile})
    get_filename_component(testname ${testsourcefile} NAME)
    #message("  New target " ${testname})
    add_executable(${testname} ${testsourcefile})
    # Make sure gstlearn is linked to each test
    target_link_libraries(${testname} gstlearn-core)
endforeach(testsourcefile ${TEST_SOURCES})


#set(TESTS_SRC Inter_Ask _AllGibbs _krige) # A completer...
#foreach(target ${TESTS_SRC})
#    add_executable(test${target} cpp/test${target}.cpp)
#    target_link_libraries(test${target} gstlearn-core)
#endforeach()

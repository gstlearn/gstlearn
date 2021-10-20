set(TESTS_SRC AllGibbs krige) # A completer...

foreach(target ${TESTS_SRC})
    add_executable(test_${target} cpp/test_${target}.cpp)
    target_link_libraries(test_${target} gstlearn-core)
endforeach()

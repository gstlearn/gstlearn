find_package(Doxygen REQUIRED)

doxygen_add_docs(doxygen ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src COMMENT "Generate man pages")
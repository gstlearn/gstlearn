find_package(Doxygen REQUIRED)

set(DOXYGEN_PROJECT_BRIEF "Geostatistics & Machine Learning toolbox")
set(DOXYGEN_MULTILINE_CPP_IS_BRIEF YES)
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_WARN_NO_PARAMDOC YES)
set(DOXYGEN_EXCLUDE ${CMAKE_SOURCE_DIR}/src/License ${CMAKE_SOURCE_DIR}/include/License ${CMAKE_SOURCE_DIR}/include/geoslib_old_f.h ${CMAKE_SOURCE_DIR}/include/geoslib_f_private.h ${CMAKE_SOURCE_DIR}/include/geoslib_d_private.h)
set(DOXYGEN_VERBATIM_HEADERS NO)
set(DOXYGEN_HTML_TIMESTAMP YES)
set(DOXYGEN_GENERATE_XML YES)
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_MACRO_EXPANSION YES)
set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
set(DOXYGEN_HIDE_UNDOC_RELATIONS NO)
set(DOXYGEN_UML_LOOK YES)
set(DOXYGEN_TEMPLATE_RELATIONS YES)

doxygen_add_docs(doxygen ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src COMMENT "Generate man pages")

# Install html and xml directories
#install(DIRECTORY ${CMAKE_BUILD_DIR}/doxygen/html/                      # HTML files from build folder
#        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/doxygen/html  # in DESTINATION/share/gstlearn/swig
#)
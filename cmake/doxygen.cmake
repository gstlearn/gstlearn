if(NOT BUILD_DOXYGEN)
  return()
endif()

# TODO : Do not regenerate doxymentation if nothing has changed in the source code
find_package(Doxygen REQUIRED)

# Configure doxyfile
set(DOXYGEN_OUTPUT_DIRECTORY doxygen)
set(DOXYGEN_PROJECT_BRIEF "Geostatistics & Machine Learning toolbox")
set(DOXYGEN_MULTILINE_CPP_IS_BRIEF YES)
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_WARN_NO_PARAMDOC YES)
set(DOXYGEN_EXCLUDE ${CMAKE_SOURCE_DIR}/include/geoslib_old_f.h
                    ${CMAKE_SOURCE_DIR}/include/geoslib_f_private.h
                    ${CMAKE_SOURCE_DIR}/include/geoslib_d_private.h)
set(DOXYGEN_VERBATIM_HEADERS NO)
set(DOXYGEN_HTML_TIMESTAMP YES)
set(DOXYGEN_GENERATE_XML YES) # For pygstlearn module
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_MAX_INITIALIZER_LINES 1000) # For very long macros
set(DOXYGEN_MACRO_EXPANSION YES)
set(DOXYGEN_EXPAND_ONLY_PREDEF NO)
set(DOXYGEN_QUIET YES)
set(DOXYGEN_HAVE_DOT NO) # Put NO to reduce generation time (keep YES for UML or better graphs)
# Uncomment if you prefer UML graphs (need DOXYGEN_HAVE_DOT YES)
#set(DOXYGEN_HIDE_UNDOC_RELATIONS NO)
#set(DOXYGEN_UML_LOOK YES)
#set(DOXYGEN_TEMPLATE_RELATIONS YES)

# Add target for generating the doxymentation
doxygen_add_docs(doxygen
                 ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src
                 COMMENT "Generate doxygen documentation")

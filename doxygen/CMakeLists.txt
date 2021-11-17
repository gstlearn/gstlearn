find_package(Doxygen REQUIRED)

# Configure doxyfile
set(DOXYGEN_PROJECT_BRIEF "Geostatistics & Machine Learning toolbox")
set(DOXYGEN_MULTILINE_CPP_IS_BRIEF YES)
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_WARN_NO_PARAMDOC YES)
set(DOXYGEN_EXCLUDE ${CMAKE_SOURCE_DIR}/src/License
                    ${CMAKE_SOURCE_DIR}/include/License
                    ${CMAKE_SOURCE_DIR}/include/geoslib_old_f.h
                    ${CMAKE_SOURCE_DIR}/include/geoslib_f_private.h
                    ${CMAKE_SOURCE_DIR}/include/geoslib_d_private.h)
set(DOXYGEN_VERBATIM_HEADERS NO)
set(DOXYGEN_HTML_TIMESTAMP YES)
set(DOXYGEN_GENERATE_XML YES) # For pygstlearn module
set(DOXYGEN_GENERATE_TREEVIEW YES)
set(DOXYGEN_MACRO_EXPANSION YES)
set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
set(DOXYGEN_HIDE_UNDOC_RELATIONS NO)
set(DOXYGEN_HAVE_DOT NO) # Put NO to reduce generation time (keep YES for UML or better graphs)
# Uncomment if you prefer UML graphs
#set(DOXYGEN_UML_LOOK YES)
#set(DOXYGEN_TEMPLATE_RELATIONS YES)

# TODO : prevent from generating twice (for build and install)

# Generate doxymentation
doxygen_add_docs(doxygen ALL
                 ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src
                 COMMENT "Generate man pages")

# Install html and xml directories
if (NOT EXCLUDE_DOXYGEN)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/                    # HTML files from build folder
          DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/doxygen/html # in DESTINATION/share/gstlearn/doxygen
  )
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/xml/                     # XML files from build folder
          DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/doxygen/xml  # in DESTINATION/share/gstlearn/doxygen
  )
endif()
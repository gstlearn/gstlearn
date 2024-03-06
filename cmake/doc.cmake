# Copy the data used by non-regression tests (in the build directory)
file(COPY doc/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc/)

# Copy the references used by demonstration and courses scripts (in the build directory)
file(COPY doc/references DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc/)

# Fix gstlearn.bib automatically
configure_file(doc/gstlearn.bib.in doc/gstlearn.bib)
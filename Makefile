# This Makefile is just a shortcut to cmake commands for Linux users only
#
# Call 'make' with one of this target:
# 
# C++ Library:
#  - shared         Build gstlearn shared library
#  - static         Build gstlearn static library
#  - build_tests    Build non-regression tests executables
#  - doxygen        Build doxygen documentation [optional]
#  - install        Install gstlearn shared library [and html doxymentation]
#  - uninstall      Uninstall gstlearn shared library [and html doxymentation]
#
# Python wrapper:
#  - python_doc     Build python package documentation [optional]
#  - python_build   Build python package [and its documentation]
#  - python_install Install python package [and its documentation]
#  - python_upload  Build python package distribution and upload to PyPi [and its documentation]
#
# No- regression tests:
#  - check_data     Execute non-regression tests (data)
#  - check_cpp      Execute non-regression tests (cpp)
#  - check_py       Execute non-regression tests (python)
#  - check          Execute non-regression tests (data + cpp + python)

# Clean:
#  - clean          Clean generated files
#  - clean_all      Clean the build directory

# You can use the following variables:
#
#  - DEBUG=1            Build the debug version of the library and tests (default =0)
#  - N_PROC=N           Use more CPUs for building procedure (default =1)
#  - BUILD_DIR=<path>   Define a specific build directory (default =build)
#  - USE_HDF5=0         To remove HDF5 support (default =1)
#

ifeq ($(DEBUG), 1)
  FLAVOR = Debug
 else
  FLAVOR = Release 
endif

ifndef BUILD_DIR
  BUILD_DIR = build
endif

ifdef N_PROC
  N_PROC_OPT = -j$(N_PROC)
endif

ifeq ($(USE_HDF5), 0)
  HDF5 = OFF
 else
  HDF5 = ON 
endif

all: shared install


.PHONY: all cmake static shared build_tests doxygen install uninstall

cmake:
	@cmake -DCMAKE_BUILD_TYPE=$(FLAVOR) -DUSE_HDF5=${HDF5} -B$(BUILD_DIR) -H.

static: cmake
	@cmake --build $(BUILD_DIR) --target static -- --no-print-directory $(N_PROC_OPT)

shared: cmake
	@cmake --build $(BUILD_DIR) --target shared -- --no-print-directory $(N_PROC_OPT)

build_tests: cmake
	@cmake --build $(BUILD_DIR) --target build_tests -- --no-print-directory $(N_PROC_OPT)

doxygen: cmake
	@cmake --build $(BUILD_DIR) --target doxygen -- --no-print-directory $(N_PROC_OPT)

install: cmake
	@cmake --build $(BUILD_DIR) --target install -- --no-print-directory $(N_PROC_OPT)

uninstall: 
	@cmake --build $(BUILD_DIR) --target uninstall -- --no-print-directory $(N_PROC_OPT)



.PHONY: python_doc python_build python_install python_upload

python_doc: cmake
	@cmake --build $(BUILD_DIR) --target python_doc -- --no-print-directory $(N_PROC_OPT)

python_build: cmake
	@cmake --build $(BUILD_DIR) --target python_build -- --no-print-directory $(N_PROC_OPT)

python_install: cmake
	@cmake --build $(BUILD_DIR) --target python_install -- --no-print-directory $(N_PROC_OPT)

python_upload: cmake
	@cmake --build $(BUILD_DIR) --target python_upload -- --no-print-directory



.PHONY: check_data check_cpp check_py check

check_data: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_data -- --no-print-directory $(N_PROC_OPT)

check_cpp: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_cpp -- --no-print-directory $(N_PROC_OPT)

check_py: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_py -- --no-print-directory $(N_PROC_OPT)

check: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check -- --no-print-directory $(N_PROC_OPT)



.PHONY: clean clean_all

clean: 
	@cmake --build $(BUILD_DIR) --target clean -- --no-print-directory $(N_PROC_OPT)

clean_all:
	@rm -rf $(BUILD_DIR)


# This Makefile is just a shortcut to cmake commands for Linux users only
#
# Call 'make' with one of this target:
# 
# C++ Library:
#  - shared         Build gstlearn shared library
#  - static         Build gstlearn static library
#  - build_tests    Build non-regression tests executables
#  - check_data     Execute non-regression tests (data)
#  - check_cpp      Execute non-regression tests (cpp)
#  - check          Execute non-regression tests (data + cpp)
#  - doxygen        Build doxygen documentation [optional]
#  - install        Install gstlearn shared library [and html doxymentation]
#  - uninstall      Uninstall gstlearn shared library [and html doxymentation]
#
# Python wrapper:
#  - python_doc     Build python package documentation [optional]
#  - python_build   Build python package [and its documentation]
#  - python_install Install python package [and its documentation]
#
# Clean:
#  - clean          Clean generated files
#  - clean_all      Clean the build directory

# You can use the following variables:
#
#  - DEBUG=1            Build the debug version of the library and tests (default =0)
#  - N_PROC=N           Use more CPUs for building procedure (default =1)
#  - BUILD_DIR=<path>   Define a specific build directory (default =build)
#

.PHONY: all cmake static shared build_tests check doxygen install uninstall

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

all: shared install

cmake:
	@cmake -DCMAKE_BUILD_TYPE=$(FLAVOR) -B$(BUILD_DIR) -H.

static: cmake
	@cmake --build $(BUILD_DIR) --target static -- --no-print-directory $(N_PROC_OPT)

shared: cmake
	@cmake --build $(BUILD_DIR) --target shared -- --no-print-directory $(N_PROC_OPT)

build_tests: cmake
	@cmake --build $(BUILD_DIR) --target build_tests -- --no-print-directory $(N_PROC_OPT)

check_data: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_data -- --no-print-directory $(N_PROC_OPT)

check_cpp: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_cpp -- --no-print-directory $(N_PROC_OPT)

check: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check -- --no-print-directory $(N_PROC_OPT)

doxygen: cmake
	@cmake --build $(BUILD_DIR) --target doxygen -- --no-print-directory $(N_PROC_OPT)

install: cmake
	@cmake --build $(BUILD_DIR) --target install -- --no-print-directory $(N_PROC_OPT)

uninstall: 
	@cmake --build $(BUILD_DIR) --target uninstall -- --no-print-directory $(N_PROC_OPT)



.PHONY: python_doc python_build python_install

python_doc: cmake
	@cmake --build $(BUILD_DIR) --target python_doc -- --no-print-directory $(N_PROC_OPT)

python_build: cmake
	@cmake --build $(BUILD_DIR) --target python_build -- --no-print-directory $(N_PROC_OPT)

python_install: cmake
	@cmake --build $(BUILD_DIR) --target python_install -- --no-print-directory $(N_PROC_OPT)



.PHONY: clean clean_all

clean: 
	@cmake --build $(BUILD_DIR) --target clean -- --no-print-directory $(N_PROC_OPT)

clean_all:
	@rm -rf $(BUILD_DIR)


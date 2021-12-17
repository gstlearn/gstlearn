# This Makefile is just a shortcut to cmake commands for Linux users only
#
# Call 'make' with one of this target:
#  - shared     Build shared library
#  - static     Build static library
#  - doxygen    Build doxygen documentation
#  - build_test Build non-regression tests executables
#  - test       Execute non-regression tests
#  - clean      Clean generated files
#  - clean_all  Clean the build directory
#  - install    Install gstlearn library
#  - uninstall  Uninstall gstlearn library
#
# You can use the following variables:
#
#  - DEBUG=1            Build the debug version of the library and tests (default =0)
#  - N_PROC=N           Use more CPUs for building procedure (default =1)
#  - BUILD_DIR=<path>   Define a specific build directory (default =build)
#

.PHONY: all cmake static shared build_test test doxygen clean clean_all

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

all: install

cmake:
	@cmake -DCMAKE_BUILD_TYPE=$(FLAVOR) -B$(BUILD_DIR) -H.

static: cmake
	@cmake --build $(BUILD_DIR) --target static -- --no-print-directory $(N_PROC_OPT)

shared: cmake
	@cmake --build $(BUILD_DIR) --target shared -- --no-print-directory $(N_PROC_OPT)

build_test: shared
	@cmake --build $(BUILD_DIR) --target build_test -- --no-print-directory $(N_PROC_OPT)

test: build_test
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target test -- --no-print-directory $(N_PROC_OPT)

doxygen: cmake
	@cmake --build $(BUILD_DIR) --target doxygen -- --no-print-directory $(N_PROC_OPT)

# TODO : I would like to install only static (but find_package(gstlearn) requires both!
# TODO : doxygen is always executed even if up-to-date
install: static shared doxygen
	@cmake --build $(BUILD_DIR) --target install -- --no-print-directory $(N_PROC_OPT)

uninstall: 
	@cmake --build $(BUILD_DIR) --target uninstall -- --no-print-directory $(N_PROC_OPT)

clean: 
	@cmake --build $(BUILD_DIR) --target clean -- --no-print-directory $(N_PROC_OPT)

clean_all:
	@rm -rf $(BUILD_DIR)


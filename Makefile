# This Makefile is just a shortcut to cmake commands for Linux users only
#
# Call 'make' with one of this target:
#  - shared     Build shared library
#  - static     Build static library
#  - doxygen    Build doxygen documentation
#  - build_test Build non-regression tests executables
#  - test       Execute non-regression tests
#  - install    Install gstlearn library
#  - uninstall  Uninstall gstlearn library
#  - clean      Clean generated files
#  - clean_all  Clean the build directory
#
# You can use the following variables:
#
#  - DEBUG=1            Build the debug version of the library and tests (default =0)
#  - N_PROC=N           Use more CPUs for building procedure (default =1)
#  - BUILD_DIR=<path>   Define a specific build directory (default =build)
#  - INSTALL_DIR=<path> Define a specific build directory (default =$HOME/gstlearn_install)
#
# Note: You must have write access to INSTALL_DIR
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

ifndef GSTLEARN_INSTALL_DIR
  # Default Linux installation folder also referenced in pygstlearn/CMakeLists.txt
  GSTLEARN_INSTALL_DIR = ${HOME}/gstlearn_install
endif

ifdef N_PROC
  N_PROC_OPT = -- -j$(N_PROC)
endif

all: static

cmake:
	@cmake -DCMAKE_BUILD_TYPE=$(FLAVOR) -DCMAKE_INSTALL_PREFIX=$(GSTLEARN_INSTALL_DIR) -B$(BUILD_DIR) -H. 

static: cmake
	@cmake --build $(BUILD_DIR) --target static $(N_PROC_OPT)

shared: cmake
	@cmake --build $(BUILD_DIR) --target shared $(N_PROC_OPT) 

build_test: shared
	@cmake --build $(BUILD_DIR) --target build_test $(N_PROC_OPT) 

test: build_test
	@cmake --build $(BUILD_DIR) --target test $(N_PROC_OPT) 

doxygen: cmake
	@cmake --build $(BUILD_DIR) --target doxygen $(N_PROC_OPT) 

# TODO : I would like to install only static (but find_package(gstlearn) requires both!
install: static shared doxygen
	@cmake --build $(BUILD_DIR) --target install

uninstall: 
	@cmake --build $(BUILD_DIR) --target uninstall

clean: 
	@cmake --build $(BUILD_DIR) --target clean $(N_PROC_OPT) 

clean_all:
	@rm -rf $(BUILD_DIR)


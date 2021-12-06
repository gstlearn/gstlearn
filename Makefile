# This Makefile is just a shortcut to cmake commands for *nix users
#
# Call 'make' with one of this target:
#  - shared     Build shared library
#  - static     Build static library
#  - doxygen    Build doxygen documentation
#  - build_test Build non-regression tests executables
#  - test       Execute non-regression tests
#  - clean      Clean generated files
#  - clean_all  Clean the build directory
#
# You can use the following variables:
#
#  - DEBUG=1          Build the debug version of the library and tests
#  - N_PROC=N         Use more CPUs for building procedure
#  - BUILD_DIR=<path> Define a specific build directory

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
  N_PROC_OPT = -- -j$(N_PROC)
endif

all: static

cmake:
	@cmake -B$(BUILD_DIR) -H. -DCMAKE_BUILD_TYPE=$(FLAVOR)

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

clean: 
	@cmake --build $(BUILD_DIR) --target clean $(N_PROC_OPT) 
	
clean_all:
	@rm -rf $(BUILD_DIR)

	


	
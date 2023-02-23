####################################################################
#
#               GSTLEARN C++ Library and Packages
#
####################################################################
#
# This Makefile is just a shortcut to cmake commands
# for make users (Linux-GCC, MacOS-clang or Windows-Rtools)
#
# Call 'make' with one of this target:
# 
# C++ Library:
#  - shared         Build gstlearn shared library
#  - static         Build gstlearn static library
#  - doxygen        Build doxygen documentation [optional]
#  - install        Install gstlearn shared library [and html doxymentation]
#  - uninstall      Uninstall gstlearn shared library [and html doxymentation]
#
# Python wrapper (only for Linux-GCC or MacOS-clang):
#  - python_doc     Build python package documentation [optional]
#  - python_build   Build python wrapper [and its documentation]
#  - python_install Install python package [and its documentation]
#  - python_upload  Build python package distribution and upload to PyPi [and its documentation]
#
# R wrapper:
#  - r_doc          Build R package documentation [optional] [TODO]
#  - r_build        Build R wrapper [and its documentation]
#  - r_install      Install R package [and its documentation]
#  - r_upload       Build R package distribution and upload to CRAN-like [and its documentation] [TODO]
#
# Non-regression tests:
#  - build_tests    Build non-regression tests C++ executables
#  - check_data     Execute non-regression tests (data)
#  - check_cpp      Execute non-regression tests (C++)
#  - check_py       Execute non-regression tests (python)
#  - check_r        Execute non-regression tests (R)
#  - check          Execute non-regression tests (data + C++ + python + R)

# Demonstration scripts:
#  - check_ipynb    Execute demonstration scripts (jupyter notebooks)
#  - check_rmd      Execute demonstration scripts (R Markdown)
#
# Clean:
#  - clean          Clean generated files
#  - clean_all      Clean the build directory
#
# You can use the following variables:
#
#  - DEBUG=1            Build the debug version of the library and tests (default =0)
#  - N_PROC=N           Use more CPUs for building procedure (default =1)
#  - BUILD_DIR=<path>   Define a specific build directory (default =build[_msys])
#  - USE_HDF5=0         To remove HDF5 support (default =1)
#
# Usage example:
#
#  make check N_PROC=2
#


ifndef USE_HDF5
  USE_HDF5 = 1
endif
ifeq ($(USE_HDF5), 1)
  USE_HDF5 = ON
 else
  USE_HDF5 = OFF 
endif

ifeq ($(OS),Windows_NT)
  # Assume MinGW (via RTools) => so MSYS Makefiles
  GENERATOR = -G"MSYS Makefiles"
else
  # Standard GNU UNIX Makefiles otherwise
  GENERATOR = -G"Unix Makefiles"
endif

ifeq ($(DEBUG), 1)
  BUILD_TYPE = Debug
 else
  BUILD_TYPE = Release 
endif

ifndef BUILD_DIR
  ifeq ($(OS),Windows_NT)
    # Assume MinGW (via RTools) => so MSYS build folder
    BUILD_DIR = build_msys
  else
    BUILD_DIR = build
  endif
endif

ifdef N_PROC
  ifeq ($(OS),Windows_NT)
    # Otherwise, tons of undefined references when compiling (don't know why)
    N_PROC_OPT = -j1
  else
    N_PROC_OPT = -j$(N_PROC)
  endif
else
  N_PROC_OPT = -j1
endif



.PHONY: all cmake cmake-python cmake-r cmake-python-r cmake-doxygen static shared build_tests doxygen install uninstall

all: shared install

cmake:
	@cmake -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DUSE_HDF5=$(USE_HDF5) -B$(BUILD_DIR) -H. $(GENERATOR)
	
cmake-python:
	@cmake -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DUSE_HDF5=$(USE_HDF5) -B$(BUILD_DIR) -H. $(GENERATOR) -DBUILD_PYTHON=ON

cmake-r:
	@cmake -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DUSE_HDF5=$(USE_HDF5) -B$(BUILD_DIR) -H. $(GENERATOR) -DBUILD_R=ON

cmake-python-r:
	@cmake -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DUSE_HDF5=$(USE_HDF5) -B$(BUILD_DIR) -H. $(GENERATOR) -DBUILD_PYTHON=ON -DBUILD_R=ON

cmake-doxygen:
	@cmake -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DUSE_HDF5=$(USE_HDF5) -B$(BUILD_DIR) -H. $(GENERATOR) -DBUILD_DOXYGEN=ON

cmake-python-doxygen:
	@cmake -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DUSE_HDF5=$(USE_HDF5) -B$(BUILD_DIR) -H. $(GENERATOR) -DBUILD_PYTHON=ON -DBUILD_DOXYGEN=ON

static: cmake
	@cmake --build $(BUILD_DIR) --target static -- --no-print-directory $(N_PROC_OPT)

shared: cmake
	@cmake --build $(BUILD_DIR) --target shared -- --no-print-directory $(N_PROC_OPT)

build_tests: cmake
	@cmake --build $(BUILD_DIR) --target build_tests -- --no-print-directory $(N_PROC_OPT)

doxygen: cmake-doxygen
	@cmake --build $(BUILD_DIR) --target doxygen -- --no-print-directory $(N_PROC_OPT)

install: cmake
	@cmake --build $(BUILD_DIR) --target install -- --no-print-directory $(N_PROC_OPT)

uninstall: 
	@cmake --build $(BUILD_DIR) --target uninstall -- --no-print-directory $(N_PROC_OPT)



.PHONY: python_doc python_build python_install python_upload

python_doc: cmake-python-doxygen
	@cmake --build $(BUILD_DIR) --target python_doc -- --no-print-directory $(N_PROC_OPT)

python_build: cmake-python
	@cmake --build $(BUILD_DIR) --target python_build -- --no-print-directory $(N_PROC_OPT)

python_install: python_build
	@cmake --build $(BUILD_DIR) --target python_install -- --no-print-directory $(N_PROC_OPT)

python_upload: python_build
	@cmake --build $(BUILD_DIR) --target python_upload -- --no-print-directory



.PHONY: r_doc r_build r_install r_upload

r_doc: cmake-r doxygen
	@echo "Target r_doc not yet implemented"

r_build: cmake-r
	@cmake --build $(BUILD_DIR) --target r_build -- --no-print-directory $(N_PROC_OPT)

r_install: r_build
	@cmake --build $(BUILD_DIR) --target r_install -- --no-print-directory $(N_PROC_OPT)

r_upload: r_build
	@echo "Target r_upload not yet implemented"



.PHONY: check_data check_cpp check_py check_r check check_ipynb check_rmd

check_data: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_data -- --no-print-directory $(N_PROC_OPT)

check_cpp: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_cpp -- --no-print-directory $(N_PROC_OPT)

check_py: cmake-python
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_py -- --no-print-directory $(N_PROC_OPT)

check_r: cmake-r
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_r -- --no-print-directory $(N_PROC_OPT)

check: cmake-python-r
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check -- --no-print-directory $(N_PROC_OPT)

check_ipynb: cmake-python
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_ipynb -- --no-print-directory $(N_PROC_OPT)

check_rmd: cmake-r
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_rmd -- --no-print-directory $(N_PROC_OPT)



.PHONY: clean clean_all

clean: 
	@cmake --build $(BUILD_DIR) --target clean -- --no-print-directory $(N_PROC_OPT)

clean_all:
	@rm -rf $(BUILD_DIR)


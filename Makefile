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
# Information:
#  - print_version  Display project name, version and date
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
#
# R wrapper:
#  - r_doc          Build R package documentation [optional]
#  - r_build        Build R wrapper [and its documentation]
#  - r_install      Install R package [and its documentation]
#
# Non-regression tests:
#  - build_tests    Build non-regression tests C++ executables
#  - check_data     Execute non-regression tests (data)
#  - check_cpp      Execute non-regression tests (C++)
#  - check_py       Execute non-regression tests (python)
#  - check_r        Execute non-regression tests (R)
#  - check          Execute non-regression tests (data + C++ + python + R)
#  - check_test_py  Execute a single python test (set $TEST variable)
#  - check_test_r   Execute a single R test      (set $TEST variable)
#  - check_test_cpp Execute a single C++ test    (set $TEST variable)
#
# Demonstration scripts:
#  - check_ipynb    Execute demonstration scripts (jupyter notebooks) and check output
#  - check_rmd      Execute demonstration scripts (R Markdown) and check output
#  - build_demos    Execute demonstration scripts (jupyter and R Markdown) and build HTMLs
#  - build_courses  Execute courses scripts (jupyter and R Markdown) and build HTMLs
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
#  - USE_HDF5=0         To remove HDF5 support (default =0)
#  - TEST=<test-target> Name of the test target to be launched (e.g. test_Model_py or test_simTub)
#  - EIGEN3_ROOT=<path> Path to Eigen3 library (optional)
#  - BOOST_ROOT=<path>  Path to Boost library (optional)
#  - LLVM_ROOT=<path>   Path to llvm compiler for MacOS only (optional)
#  - SWIG_EXEC=<path>   Path to swig executable (optional)
#
# Usage example:
#
#  make check N_PROC=2
#

ifdef USE_HDF5
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
  # Set OS also for Linux or Darwin
  OS := $(shell uname -s)
endif

ifeq ($(OS),Darwin)
  ifndef LLVM_ROOT
  	LLVM_ROOT = /opt/homebrew
  endif
  # Particular clang compiler for supporting OpenMP
  CC_CXX = CC=$(LLVM_ROOT)/opt/llvm/bin/clang CXX=$(LLVM_ROOT)/opt/llvm/bin/clang++
else
  CC_CXX = 
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

CMAKE_DEFINES := -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DUSE_HDF5=$(USE_HDF5)
ifdef SWIG_EXEC
  CMAKE_DEFINES := $(CMAKE_DEFINES) -DSWIG_EXECUTABLE=$(SWIG_EXEC)
endif
ifdef EIGEN3_ROOT
  CMAKE_DEFINES := $(CMAKE_DEFINES) -DEigen3_ROOT=$(EIGEN3_ROOT)
endif
ifdef BOOST_ROOT
  CMAKE_DEFINES := $(CMAKE_DEFINES) -DBoost_ROOT=$(BOOST_ROOT)
endif

.PHONY: all cmake cmake-python cmake-r cmake-python-r cmake-doxygen print_version static shared build_tests doxygen install uninstall

all: shared install

cmake:
	@$(CC_CXX) cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES)

cmake-python:
	@$(CC_CXX) cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_PYTHON=ON

cmake-r:
	@$(CC_CXX) cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_R=ON

cmake-python-r:
	@$(CC_CXX) cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_PYTHON=ON -DBUILD_R=ON

cmake-doxygen:
	@$(CC_CXX) cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_DOXYGEN=ON

cmake-python-doxygen:
	@$(CC_CXX) cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_PYTHON=ON -DBUILD_DOXYGEN=ON

cmake-r-doxygen:
	@$(CC_CXX) cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_R=ON -DBUILD_DOXYGEN=ON

print_version: cmake
	@cmake --build $(BUILD_DIR) --target print_version -- --no-print-directory

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



.PHONY: python_doc python_build python_install

python_doc: cmake-python-doxygen
	@cmake --build $(BUILD_DIR) --target python_doc -- --no-print-directory $(N_PROC_OPT)

python_build: cmake-python
	@cmake --build $(BUILD_DIR) --target python_build -- --no-print-directory $(N_PROC_OPT)

python_install: cmake-python
	@cmake --build $(BUILD_DIR) --target python_install -- --no-print-directory $(N_PROC_OPT)


.PHONY: r_doc r_build r_install

r_doc: cmake-r #cmake-r-doxygen
	@cmake --build $(BUILD_DIR) --target r_doc -- --no-print-directory $(N_PROC_OPT)

r_build: cmake-r #cmake-r-doxygen
	@cmake --build $(BUILD_DIR) --target r_build -- --no-print-directory $(N_PROC_OPT)

r_install: cmake-r #cmake-r-doxygen
	@cmake --build $(BUILD_DIR) --target r_install -- --no-print-directory $(N_PROC_OPT)


.PHONY: check_data check_cpp check_py check_r check check_ipynb check_rmd check_test_cpp check_test_py check_test_r build_demos build_courses

check_data: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_data -- --no-print-directory $(N_PROC_OPT)

check_cpp: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_cpp -- --no-print-directory $(N_PROC_OPT)

check_py: cmake-python
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_py -- --no-print-directory $(N_PROC_OPT)

check_r: cmake-r #cmake-r-doxygen
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_r -- --no-print-directory $(N_PROC_OPT)

check: cmake-python-r
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check -- --no-print-directory $(N_PROC_OPT)

check_ipynb: cmake-python
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_ipynb -- --no-print-directory $(N_PROC_OPT)

check_rmd: cmake-r #cmake-r-doxygen
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_rmd -- --no-print-directory $(N_PROC_OPT)

check_test_cpp: cmake
	@cd $(BUILD_DIR); make $(TEST); CTEST_OUTPUT_ON_FAILURE=1 ctest -R $(TEST)

check_test_py: cmake-python
	@cd $(BUILD_DIR); make prepare_check_py; CTEST_OUTPUT_ON_FAILURE=1 ctest -R $(TEST)

check_test_r: cmake-r-doxygen
	@cd $(BUILD_DIR); make prepare_check_r; CTEST_OUTPUT_ON_FAILURE=1 ctest -R $(TEST)

dump_test_cpp: cmake
	@cd $(BUILD_DIR); make $(TEST); "tests/cpp/$(BUILD_TYPE)/$(TEST)" dummy

build_demos: cmake-python-r
	@cmake --build $(BUILD_DIR) --target build_demos -- --no-print-directory $(N_PROC_OPT)

build_courses: cmake-python-r
	@cmake --build $(BUILD_DIR) --target build_courses -- --no-print-directory $(N_PROC_OPT)


.PHONY: clean clean_all

clean: 
	@cmake --build $(BUILD_DIR) --target clean -- --no-print-directory $(N_PROC_OPT)

clean_all:
	@rm -rf $(BUILD_DIR)


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
#  - BUILD_PYTHON=1     Configure cmake to build python wrapper (default =0, see target python_*)
#  - BUILD_R=1          Configure cmake to build R wrapper (default =0, see target r_*)
#  - BUILD_DOC=1        Configure cmake to build documentation (default =0)
#  - TEST=<test-target> Name of the test target to be launched (e.g. test_Model_py or test_simTub)
#  - ASAN=1             Build with Address Sanitizer (default =0)
#  - USE_HDF5=0         To remove HDF5 support (default =1)
#  - NO_INTERNET=1      To prevent python pip from looking for dependencies through Internet
#                       (useful when there is no Internet available) (default =0)
#  - NO_BYTE_COMPILE=1  To prevent R from byte-compiling the R package during the installation
#                       (useful to accelerate installation process) (default =0) 
#  - EIGEN3_ROOT=<path> Path to Eigen3 library (optional)
#  - BOOST_ROOT=<path>  Path to Boost library (optional)
#  - NLOPT_ROOT=<path>  Path to NLopt library (optional)
#  - LLVM_ROOT=<path>   Path to llvm compiler for MacOS only (optional)
#  - SWIG_EXEC=<path>   Path to swig executable (optional)
#
# Usage example:
#
#  make check N_PROC=2
#

ifeq ($(BUILD_DOC), 1)
  BUILD_DOC = ON
 else
  BUILD_DOC = OFF 
endif

ifeq ($(NO_INTERNET), 1)
  NO_INTERNET = ON
 else
  NO_INTERNET = OFF 
endif

ifeq ($(NO_BYTE_COMPILE), 1)
  NO_BYTE_COMPILE = ON
 else
  NO_BYTE_COMPILE = OFF
endif

ifeq ($(USE_HDF5), 0)
  USE_HDF5 = OFF
 else
  USE_HDF5 = ON
endif

ifeq ($(OS),Windows_NT)
  # Assume MinGW (via RTools) => so MSYS Makefiles
  GENERATOR = -G"MSYS Makefiles"
else
  ifeq (, $(shell which ninja))
    # Standard GNU UNIX Makefiles otherwise
    GENERATOR = -G "Unix Makefiles"
  else
    # Standard GNU UNIX Makefiles otherwise
    GENERATOR = -G "Ninja"
  endif
  # Set OS also for Linux or Darwin
  OS := $(shell uname -s)
endif

ifeq ($(DEBUG), 1)
  BUILD_TYPE = Debug
 else
  BUILD_TYPE = Release
endif

ifeq ($(ASAN), 1)
  BUILD_ASAN = ON
 else
  BUILD_ASAN = OFF
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
    N_PROC_OPT = -j1 | tee /dev/null
  else
    N_PROC_OPT = -j$(N_PROC) | tee /dev/null
  endif
else
  N_PROC_OPT = -j1 | tee /dev/null
endif
# Add  "| tee /dev/null" because Ninja prints output in a single line :
# https://stackoverflow.com/questions/46970462/how-to-enable-multiline-logs-instead-of-single-line-progress-logs

CMAKE_DEFINES := -DCMAKE_BUILD_TYPE=$(BUILD_TYPE) -DUSE_HDF5=$(USE_HDF5) -DBUILD_ASAN=$(BUILD_ASAN) -DBUILD_TESTING=ON
ifdef SWIG_EXEC
  CMAKE_DEFINES := $(CMAKE_DEFINES) -DSWIG_EXECUTABLE=$(SWIG_EXEC)
endif
ifdef EIGEN3_ROOT
  CMAKE_DEFINES := $(CMAKE_DEFINES) -DEigen3_ROOT=$(EIGEN3_ROOT)
endif
ifdef BOOST_ROOT
  CMAKE_DEFINES := $(CMAKE_DEFINES) -DBoost_ROOT=$(BOOST_ROOT)
endif
ifdef NLOPT_ROOT
  CMAKE_DEFINES := $(CMAKE_DEFINES) -DNLopt_ROOT=$(NLOPT_ROOT)
endif

.PHONY: all cmake cmake-python cmake-r cmake-python-r cmake-doxygen print_version static shared build_tests doxygen install uninstall

all: shared install

cmake:
	@cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_DOC=$(BUILD_DOC)

cmake-python:
	@cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_DOC=$(BUILD_DOC) -DBUILD_PYTHON=ON -DNO_INTERNET=$(NO_INTERNET)

cmake-r:
	@cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_DOC=$(BUILD_DOC) -DBUILD_R=ON -DNO_BYTE_COMPILE=$(NO_BYTE_COMPILE)

cmake-python-r:
	@cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_DOC=$(BUILD_DOC) -DBUILD_PYTHON=ON -DBUILD_R=ON -DNO_INTERNET=$(NO_INTERNET) -DNO_BYTE_COMPILE=$(NO_BYTE_COMPILE)

cmake-doxygen:
	@cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_DOC=ON

cmake-python-doxygen:
	@cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_DOC=ON -DBUILD_PYTHON=ON -DNO_INTERNET=$(NO_INTERNET)

cmake-r-doxygen:
	@cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_DOC=ON -DBUILD_R=ON -DNO_BYTE_COMPILE=$(NO_BYTE_COMPILE)

cmake-python-r-doxygen:
	@cmake -B$(BUILD_DIR) -S. $(GENERATOR) $(CMAKE_DEFINES) -DBUILD_DOC=ON -DBUILD_PYTHON=ON -DBUILD_R=ON -DNO_INTERNET=$(NO_INTERNET) -DNO_BYTE_COMPILE=$(NO_BYTE_COMPILE)

print_version: cmake
	@cmake --build $(BUILD_DIR) --target print_version --

static shared build_tests doxygen install uninstall: cmake-doxygen
	@cmake --build $(BUILD_DIR) --target $@ -- $(N_PROC_OPT)


.PHONY: python_doc python_build python_install

python_doc: cmake-python-doxygen
	@cmake --build $(BUILD_DIR) --target python_doc -- $(N_PROC_OPT)

python_build: cmake-python
	@cmake --build $(BUILD_DIR) --target python_build -- $(N_PROC_OPT)

python_install: cmake-python
	@cmake --build $(BUILD_DIR) --target python_install -- $(N_PROC_OPT)


.PHONY: r_doc r_build r_install

r_doc: cmake-r-doxygen
	@cmake --build $(BUILD_DIR) --target r_doc -- $(N_PROC_OPT)

r_build: cmake-r
	@cmake --build $(BUILD_DIR) --target r_build -- $(N_PROC_OPT)

r_install: cmake-r
	@cmake --build $(BUILD_DIR) --target r_install -- $(N_PROC_OPT)


.PHONY: check_data check_cpp check_py check_r check check_ipynb check_rmd check_test_cpp check_test_py check_test_r build_demos build_courses

check_data: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_data -- $(N_PROC_OPT)

check_cpp: cmake
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_cpp -- $(N_PROC_OPT)

check_py: cmake-python
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_py -- $(N_PROC_OPT)

check_r: cmake-r
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_r -- $(N_PROC_OPT)

check: cmake-python-r
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check -- $(N_PROC_OPT)

check_ipynb: cmake-python
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_ipynb -- $(N_PROC_OPT)

check_rmd: cmake-r
	@CTEST_OUTPUT_ON_FAILURE=1 cmake --build $(BUILD_DIR) --target check_rmd -- $(N_PROC_OPT)

check_test_cpp: cmake
	@cd $(BUILD_DIR); make $(TEST); CTEST_OUTPUT_ON_FAILURE=1 ctest -R '^$(TEST)$$'

check_test_py: cmake-python
	@cd $(BUILD_DIR); make python_install; make prepare_check_py; make prepare_check_ipynb; CTEST_OUTPUT_ON_FAILURE=1 ctest -R '^$(TEST)$$'

check_test_r: cmake-r
	@cd $(BUILD_DIR); make r_install; make prepare_check_r; make prepare_check_rmd; CTEST_OUTPUT_ON_FAILURE=1 ctest -R '^$(TEST)$$'

dump_test_cpp: cmake
	@cd $(BUILD_DIR); make $(TEST); "tests/cpp/$(BUILD_TYPE)/$(TEST)" dummy

build_demos: cmake-python-r
	@cmake --build $(BUILD_DIR) --target build_demos -- $(N_PROC_OPT)

build_courses: cmake-python-r
	@cmake --build $(BUILD_DIR) --target build_courses -- $(N_PROC_OPT)

.PHONY: scan_build clang_tidy clang_check

scan_build: clean_all
	@scan-build cmake -B$(BUILD_DIR) -S. $(CMAKE_DEFINES) -DBUILD_TESTING=ON
	@CCACHE_DISABLE=1 scan-build -V --exclude Eigen --exclude 3rd-party cmake --build $(BUILD_DIR)

clang_tidy: clean_all
	@cmake -B$(BUILD_DIR) -S. $(CMAKE_DEFINES) -DBUILD_TESTING=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
	@run-clang-tidy -p $(BUILD_DIR)

clang_check: clean_all
	@cmake -B$(BUILD_DIR) -S. $(CMAKE_DEFINES) -DBUILD_TESTING=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
	@git ls-files | grep -E "src.*\.cpp|include.*\.hpp" | xargs -n1 -P$(shell nproc) clang-check -p $(BUILD_DIR)

.PHONY: clean clean_all

clean: 
	@cmake --build $(BUILD_DIR) --target clean -- $(N_PROC_OPT)

clean_all:
	@rm -rf $(BUILD_DIR)


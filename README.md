## Overview

*gstlearn* is the new cross-platform Geostatistics C++ Library proposed by MINES Paris - PSL University. It offers to C++ users **all famous Geostatistical methodologies** developed and/or invented by the Geostatistic Team of the [Geosciences Research Center](https://www.geosciences.minesparis.psl.eu/).<br/>
Copyright (c) MINES Paris / PSL University

The name 'gstlearn' stands for several purposes:

* GeoSTatistics & Machine Learning Library
* Geostatistical Spatio-Temporal Learning
* Learning Geosciences & Spatio-Temporal Models

*gstlearn* comes in different forms:

* A C++ Library (this repository): [https://github.com/gstlearn/gstlearn](https://github.com/gstlearn/gstlearn)
* A Python Package: [https://github.com/gstlearn/gstlearn/tree/main/python](https://github.com/gstlearn/gstlearn/tree/main/python)
* A R Package: TODO: coming soon (meanwhile, you may use [RGeostats R package](http://cg.ensmp.fr/rgeostats))

## References

The *gstlearn* C++ Library is the direct successor of the Geoslib C/C++ Library which was proposed through the [RGeostats R package](http://cg.ensmp.fr/rgeostats).

The *gstlearn* C++ Library is developed by the [Geostatistics Group](https://www.geosciences.minesparis.psl.eu/en/presentation/geostatistics) of the [Geosciences Center](https://www.geosciences.minesparis.psl.eu) ([MINES PariTech](https://mines-paristech.eu/) - [PSL University](https://psl.eu/en) - France)

When using the *gstlearn* C++ Library, please use the citation from [doc/gstlearn.bib](doc/gstlearn.bib)

## Requirements

This library has been successfully tested with Ubuntu 16/18/20 LTS and Windows 10 (MacOS: not tested).
For compiling and installing *gstlearn* C++ Library, the following tools must be available (See [required tools installation](#required-tools-installation) instructions below):

* [Git](https://git-scm.com/downloads) client
* [CMake](https://cmake.org/download) tool 3.19 or higher
* A C++ compiler:
  * Linux/MacOS:
   * [GCC](https://gcc.gnu.org) compiler 5.4 or higher
  * Windows:
   * [Microsoft Visual Studio C++](https://visualstudio.microsoft.com/fr/vs/features/cplusplus/) 14 or higher
   * [MinGW](https://wiki.python.org/moin/WindowsCompilers#GCC_-_MinGW-w64_.28x86.2C_x64.29) 7 or higher
* [Doxygen](https://www.doxygen.nl/download.html) 1.8.3 or higher
* [Boost](https://www.boost.org/users/download) header files
* [HDF5](https://www.hdfgroup.org/solutions/hdf5/) C & C++ library and header files 1.8 or higher
  
## Get the sources

For getting the sources files, just clone the github repository:

    git clone https://github.com/gstlearn/gstlearn.git
    cd gstlearn

Notes:

* In the following, all instructions must be executed from a command prompt inside this *root* directory (thus the last command `cd gstlearn`)

## Library compilation & installation

For compiling and installing the *gstlearn* C++ shared Library, execute the following instructions. Please note that you can choose another destination folder (currently named *build*).

### GCC, Clang, MinGW, ...

    cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Release
    cmake --build build --target shared
    cmake --build build --target install

or for those who prefer a single command line

    mkdir -p build & cd build & cmake .. & make shared & make install

or even faster

    make

Note:

* Using MingGW on a Windows where Visual Studio is also installed may need to add `-G "MSYS Makefiles"` in the first command.

### Microsoft Visual Studio, XCode, ...

    cmake -Bbuild -H.
    cmake --build build --target shared --config Release
    cmake --build build --target install --config Release

Note:

* Using Visual Studio on a Windows where MingGW is also installed may need to add `-G "Visual Studio 16 2019"` in the first command (adapt version).
  
### Important Notes

Notes:

* The default installation directory named *gstlearn_install* is located in your *Home*. If you want to change it, you can either:
    * Define the `GSTLEARN_INSTALL_DIR` environment variable or
    * Add `-DGSTLEARN_INSTALL_DIR=<path/of/gstlearn/install/dir>` to the first cmake command above
* If you want to build and install the *Debug* version, you must replace `Release` by `Debug` above
* The *static* version of the library is mandatory for creating [python package](https://github.com/gstlearn/gstlearn/tree/main/python)
* Only the *shared* library (built by default) is installed.
* You may need to precise the location of Boost or HDF5 installation directory (which contain *include* and *lib* folders). In that case, add the following in the first command above:
    * `-DBoost_ROOT=<path/to/boost>`
    * `-DHDF5_ROOT=<path/to/hdf5>`

## Usage

Please, look at *tests* C++ code in order to learn how to use the gstlearn C++ library.

## Required tools installation
### Linux (Ubuntu):

    sudo apt install git
    sudo apt install cmake
    sudo apt install doxygen
    sudo apt install libboost-all-dev
    sudo apt install libhdf5-dev

Notes:

* Under Linux, the GCC compiler and GNU make is already installed
* If your Linux distribution repositories don't provide minimum required versions, please install the tools manually (see provider website)

### MacOS:

    brew install git
    brew install cmake
    brew install doxygen
    brew install libboost-all-dev
    brew install libhdf5-dev

Notes:

* This is currently not tested - above packages may not exist
* Under MacOS, the GCC (or Clang) compiler and GNU make is already installed
* If your MacOS repositories don't provide minimum required versions, please install the tools manually (see provider website)
  
### Windows:

Download and install the following tools:

* Git client [from here](https://gitforwindows.org) (Use default options during installation)
* CMake tool [from here](https://cmake.org/download) (Check the 'Add CMake to the Path' option during installation)
* Doxygen 1.8.3+ [from here](https://www.doxygen.nl/download.html) (Install in the directory *C:\\doxygen* for example)

Notes:

* You must restart your computer after installing these requirements
* The *Path* environment variable must be updated to make *doxygen.exe* available (follow [this guide](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho) to add *C:\\doxygen\\bin* folder in the *Path* variable and restart Windows)

#### Microsoft Visual Studio

Download and install the following tools:

* Microsoft Visual Studio C++ 14+ [from here](https://visualstudio.microsoft.com/fr/vs/features/cplusplus/)
* Boost library [from here](https://www.boost.org/users/download) (Download and extract the zip file in *C:\\local\\* directory for example)
* HDF5 library [from here](https://www.hdfgroup.org/downloads/hdf5) (Download the pre-built binaries (zip), extract the zip file and execute the installer using default options)

#### MingGW

Download and install the following tools:

* Rtools4 [from here](https://cran.r-project.org/bin/windows/Rtools/rtools40.html)
  
Notes:

* You must restart your computer after installing these requirements
* Rtools is not the unique way to install MinGW, but it is the best way to create R packages on Windows.

Then, from a Windows command prompt, execute following instructions:

    pacman -S mingw-w64-x86_64-hdf5
    pacman -S mingw-w64-x86_64-boost

## Development
### Non-regression tests
#### GCC, Clang, MinGW, ...

To launch non-regression tests, execute the following command:

    cmake --build build --target build_tests
    cmake --build build --target check

or for those who prefer a single command line

    make check

#### Microsoft Visual Studio, XCode, ...

To build and launch non-regression tests, execute the following commands:

    cmake --build build --target build_tests --config Release
    cmake --build build --target check --config Release

Notes:

* If you want to run the *Debug* version of the tests, you must replace `Release` by `Debug` above
* The *check* target brings some required runtime customization, so do not use the standard *ctest* command
  
### Clean

To clean (partially) the build, execute the following command:

    cmake --build build --target clean

Notes:

* If you really want to clean all files generated by CMake, you can remove *build* directory content by hand. Linux or MinGW users may use the clean_all target from the source directory.

### Uninstall

To uninstall all the installed files (only the files, not the directories), execute this command:

    cmake --build build --target uninstall

### Generate the documentation

The doxygen HTML documentation is optional (not included in the installation by default). If you want to generate it, execute the command:

    cmake --build build --target doxygen

The documentation is then available by opening the following HTML file with your favorite web-browser:

    firefox build/doxygen/html/index.html

---

## License
MIT
2022 Team gstlearn

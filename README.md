## Overview

**gstlearn** is the new cross-platform Geostatistics C++ library proposed by MINES Paris - PSL University. It offers to C++ users **all famous Geostatistical methodologies** developed and/or invented by the Geostatistic Team of the [Geosciences Research Center](https://www.geosciences.minesparis.psl.eu/).<br/>
Copyright (c) MINES PARIS / PSL University

The name 'gstlearn' stands for several purposes:

* GeoSTatistics & Machine Learning Library
* Geostatistical Spatio-Temporal Learning
* Learning Geosciences & Spatio-Temporal Models

*gstlearn* comes in various forms:

* A C++ library: [https://github.com/gstlearn/gstlearn](https://github.com/gstlearn/gstlearn)
* A Python package: [https://github.com/gstlearn/gstlearn/tree/main/python](https://github.com/gstlearn/gstlearn/tree/main/python)
* A R package: [https://github.com/gstlearn/gstlearn/tree/main/r](https://github.com/gstlearn/gstlearn/tree/main/r)

If you simply want to install the Python or R package for gstlearn you should look at the corresponding sub folders.

## References

The *gstlearn* C++ library is the direct successor of the Geoslib C/C++ library which was proposed through the [RGeostats R package](http://cg.ensmp.fr/rgeostats).

The *gstlearn* C++ library is developed by the [Geostatistics Group](https://www.geosciences.minesparis.psl.eu/en/presentation/geostatistics) of the [Geosciences Center](https://www.geosciences.minesparis.psl.eu) ([MINES Paris](https://mines-paristech.eu/) - [PSL University](https://psl.eu/en) - France)

When using the *gstlearn* C++ library, please use the citation from [doc/gstlearn.bib](doc/gstlearn.bib.in)

The *gstlearn* C++ library is a derivative work based on the *swigex* project: [https://github.com/fabien-ors/swigex](https://github.com/fabien-ors/swigex)

## Requirements

This library has been successfully tested with Ubuntu 16/18/20 LTS and Windows 10 (MacOS: not tested).
For compiling and installing *gstlearn* C++ library, the following tools must be available (See [required tools installation](#required-tools-installation) instructions below):

* Git client
* CMake tool 3.20 or higher
* A C++ compiler among:
  * Linux/MacOS:
    * GCC 5.4 or higher
    * Clang 12 or higher (not tested)
  * Windows:
    * Python users: Microsoft Visual Studio C++ 14 or higher
    * R users: MinGW 7 (RTools) or higher
* Boost header files
* Doxygen [Optional] 1.8.3 or higher
* HDF5 [Optional] C++ library and header files 1.8 or higher

See [required tools installation](#required-tools-installation) instructions below

## Get the sources

For getting the sources files, just clone the github repository:

    git clone https://github.com/gstlearn/gstlearn.git
    cd gstlearn

Notes:

* In the following, all instructions must be executed from a command prompt inside this *root* directory (thus the last command `cd gstlearn` above)

## C++ Library Compilation & Installation

For compiling and installing the *gstlearn* C++ shared library, execute the following instructions. Please note that you can choose another destination folder (currently named *build*).

If you only want to use Python or R packages, you should switch to corresponding README:
* Python: [https://github.com/gstlearn/gstlearn/tree/main/python](https://github.com/gstlearn/gstlearn/tree/main/python)
* R: [https://github.com/gstlearn/gstlearn/tree/main/r](https://github.com/gstlearn/gstlearn/tree/main/r)

### GCC, Clang, MinGW, ...

...or any other single configuration compiler:

    cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Release
    cmake --build build --target shared
    cmake --build build --target install

or for those who prefer a single command line:

    mkdir -p build & cd build & cmake .. & make shared & make install

or even faster:

    make

Note:

* Using MinGW on a Windows where another compiler is also installed may need to add `-G "MSYS Makefiles"` in the first cmake command above.

### Microsoft Visual Studio, ...

...or any other multiple configurations compiler:

    cmake -Bbuild -H.
    cmake --build build --target shared --config Release
    cmake --build build --target install --config Release

Note:

* Using Visual Studio on a Windows where another compiler is also installed may need to add `-G "Visual Studio 16 2019"` in the first command (adapt version).
  
### Important Notes

Notes:

* Currently, **HDF5 is not supported** when compiling *gstlearn* C++ library **under Windows**. *gstlearn* won't link against HDF5 and GibbsMMulti::setFlagStoreInternal(false) feature won't be available.
* The default installation directory named *gstlearn_install* is located in your *Home*. If you want to change it, you can either:
    * Define the `GSTLEARN_INSTALL_DIR` environment variable or
    * Add `-DGSTLEARN_INSTALL_DIR=<path/of/gstlearn/install/dir>` to the first cmake command above
* If you want to build and install the *Debug* version, you must replace `Release` by `Debug` above
* If you don't want HDF5 support, add `-DUSE_HDF5=OFF` to the first cmake command above. If you use the shortcut Makefile, you can use `USE_HDF5=0` after the `make` command
* Only the *shared* library (built by default) is installed when compiling *gstlearn* C++ library. If you want to compile the *static* version, you must replace *shared* by *static* target above.
* You may need to precise the location of Boost, HDF5 or Doxygen installation directory. In that case, add the following variables in the first cmake command above:
    * `-DBoost_ROOT=<path/to/boost>`
    * `-DHDF5_ROOT=<path/to/hdf5>`
    * `-DDoxygen_ROOT=<path/to/doxygen>`

## Usage

Please, look at *tests* [C++ code](https://github.com/gstlearn/gstlearn/tree/main/tests) in order to learn how to use the *gstlearn* C++ library. You can generate the source code documentation using [Doxygen](#generate-the-documentation).

## Required Tools Installation

These tools are needed for compiling the *gstlearn* C++ library. Please note that HDF5 and Doxygen installation are optional.

### Linux (Ubuntu)

    sudo apt install git
    sudo apt install cmake
    sudo apt install doxygen
    sudo apt install libboost-all-dev
    sudo apt install libhdf5-dev

Notes:

* Under Linux, the GCC compiler and GNU make is already installed
* If your Linux distribution repositories don't provide minimum required versions, please install the tools manually (see provider website)

### MacOS

    brew install git
    brew install cmake
    brew install doxygen
    brew install libboost-all-dev
    brew install libhdf5-dev

Notes:

* These instructions for MacOS are currently not tested - above packages may not exist
* Under MacOS, the GCC (or Clang) compiler and GNU make is already installed
* If your MacOS repositories don't provide minimum required versions, please install the tools manually (see provider website)
  
### Windows - Microsoft Visual Studio

These requirements are also recommended to people who wants to compile *gstlearn* Python package. If you want to compile *gstlearn* R package under Windows, you should look at the next section.

Download and install the following tools using default options during installation:

* Git client [from here](https://gitforwindows.org) (*Setup program* [exe])
* CMake tool [from here](https://cmake.org/download) (*Windows Installer* [msi], check the *'Add CMake to the system PATH for all users'* option during installation)
* Microsoft Visual Studio (Community) [from here](https://visualstudio.microsoft.com/fr/thank-you-downloading-visual-studio/?sku=Community&channel=Release&version=VS2022&source=VSLandingPage&cid=2030&passive=false) (*VisualStudioSetup.exe* - only select the *Visual Studio Desktop C++* component)
* Boost library [from here](https://www.boost.org/users/download) (*Archive file* [zip] to be extracted in a folder of your own - and remind that folder)
* HDF5 library (optional) [from here](https://www.hdfgroup.org/downloads/hdf5) (*Pre-built binaries* [zip] to be extracted, then, execute the *installer* [msi] - and remind the installation folder, we assume it is `C:\Program Files\HDF_Group\HDF5\1.12.2`)
* Doxygen (optional) [from here](https://www.doxygen.nl/download.html) (*Binary distribution* [setup.exe] - remind the installation folder, we assume it is `C:\Program Files\doxygen`)

Notes:

* The *Path* environment variable (*System variables*) must be updated to make *doxygen.exe* available in the batch command line (follow [this guide](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10) to add *bin* directory from the *Doxygen* installation folder in the *Path* variable and restart Windows)
* You must restart your computer after installing these requirements

### Windows - MinGW (via RTools):

These requirements are also recommended to people who wants to compile *gstlearn* R package. If you want to compile *gstlearn* Python package under Windows, you should look at the previous section. This is not the only way to install MinGW. But using MinGW provided with RTools permits us to also handle *gstlearn* R package compilation.

#### Install R and RTools

Download and install the following tools using default option during installation:

* R [from here](https://cran.r-project.org/bin/windows/base) (*Setup program* [exe] - remind the installation folder, assume it is `C:\Program Files\R\R-4.2.2`)
* RTools [from here](https://cran.r-project.org/bin/windows/Rtools) (RTools *Installer* [exe] - remind the installation folder, assume it is `C:\rtools42`)

Notes:

* Choose the corresponding RTools version according to the R version installed
* Instructions in this section are **valid since R v4.2** (for older versions please contact us)
* RTools is not the unique way to install MinGW on Windows, but it is our preferred way as we can handle R packages compilation
* The *Path* environment variable (*System variables*) must be updated to make *R.exe* available in the batch command line (follow [this guide](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10) to add `bin` directory from the *R* installation folder in the *Path* variable and restart Windows, ie: `C:\Program Files\R\R-4.2.2\bin`)

#### Add MSYS2 Required Packages

1. Edit the `etc/pacman.conf` file in the RTools installation directory (ie: `C:\rtools42`) by changing the `SigLevel` variable to `Never` (otherwise, *git* cannot be installed using *pacman*):

````
SigLevel=Never
````

2. Edit the `mingw64.ini` file in the RTools installation directory (ie: `C:\rtools42`) by un-commenting the following line (remove '#' character at the beginning):

````
MSYS2_PATH_TYPE=inherit
````

3. Launch *mingw64.exe* in RTools installation directory (ie: `C:\rtools42`) and pin the icon to the task bar

4. From the *mingw64* shell command prompt, execute following instructions (HDF5 and Doxygen are optional):

````
pacman -Sy git
pacman -Sy mingw-w64-x86_64-cmake
pacman -Sy mingw-w64-x86_64-gcc
pacman -Sy mingw-w64-x86_64-boost
pacman -Sy mingw-w64-x86_64-hdf5
pacman -Sy mingw-w64-x86_64-doxygen
````

## Development

### Execute Non-regression Tests

The `check.*` targets brings some required runtime customization, so do not use the standard *ctest* command for triggering the non-regression tests.

To build and launch non-regression tests, execute the following commands:

#### GCC, Clang, MinGW, ...

...or any other single configuration compiler:

    cmake --build build --target build_tests
    cmake --build build --target check_cpp
    cmake --build build --target check_data
    
or even faster:

    make check_cpp
    make check_data

#### Microsoft Visual Studio, ...

...or any other multiple configurations compiler:

    cmake --build build --target build_tests --config Release
    cmake --build build --target check_cpp --config Release
    cmake --build build --target check_data --config Release
  
### Clean

To clean (partially) the build, execute the following command:

    cmake --build build --target clean

Notes:

* If you really want to clean all files generated by CMake, you can remove *build* directory content by hand. Linux, MacOS or MinGW users may use the `clean_all` target from the shortcut Makefile inside the top level directory:

````
make clean_all
````

### Uninstall the Library

To uninstall all the installed files (only the files, not the directories), execute this command:

    cmake --build build --target uninstall

or faster:

    make uninstall

### Generate the Documentation

The Doxygen HTML documentation is optional (not included in the installation by default). If you want to generate it, execute the command:

    cmake --build build --target doxygen

or faster:

    make doxygen

The documentation is then available by opening the following HTML file with your favorite web-browser:

    firefox build/doxygen/html/index.html

---

## License
BSD
2022 Team gstlearn

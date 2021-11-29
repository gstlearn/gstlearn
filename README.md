## Overview
*gstlearn* is the new cross-platform Geostatistics C++ Library proposed by MINES Paristech - PSL University. It offers to C++ users **all famous Geostatistical methodologies** developed and/or invented by the Geostatistic Team of the [Geosciences Research Center](https://www.geosciences.minesparis.psl.eu/)!<br/>
Copyright (c) MINES Paristech / PSL University

The name 'gstlearn' stands for several purposes:
  * GeoSTatistics & Machine Learning Library
  * Geostatistical Spatio-Temporal Learning
  * Learning Geosciences & Spatio-Temporal Models

*gstlearn* comes in different forms:
  * A C++ Library (this repository): [https://github.com/gstlearn/gstlearn](https://github.com/gstlearn/gstlearn)
  * A Python Package: [https://github.com/gstlearn/pygstlearn](https://github.com/gstlearn/pygstlearn)
  * A R Package: TODO: coming soon (meanwhile, you may use [RGeostats R package](http://cg.ensmp.fr/rgeostats))

## References
The *gstlearn* C++ Library is the direct successor of the Geoslib C/C++ Library which was proposed through the [RGeostats R package](http://cg.ensmp.fr/rgeostats).

The *gstlearn* C++ Library is developed by the [Geostatistics Group](https://www.geosciences.minesparis.psl.eu/en/presentation/geostatistics) of the [Geosciences Center](https://www.geosciences.minesparis.psl.eu) ([MINES PariTech](https://mines-paristech.eu/) - [PSL University](https://psl.eu/en) - France)

When using the *gstlearn* C++ Library, please use the citation from [doc/gstlearn.bib](doc/gstlearn.bib)

## Requirements
This package has been successfully tested with Ubuntu 16.04 LTS, Ubuntu 18.04 LTS and Windows 10 (MacOS: not tested).
For compiling and installing *gstlearn* C++ Library, the following tools must be available (See [required tools installation](#required-tools-installation) instructions below):
  * [Git](https://git-scm.com/downloads) client
  * [CMake](https://cmake.org/download) tool 3.19 or higher
  * A C++ compiler:
    * Linux/MacOS:
      * [GCC](https://gcc.gnu.org) compiler 5.4 or higher
    * Windows:
      * [Microsoft Visual C++ Compiler](https://visualstudio.microsoft.com/visual-cpp-build-tools) 14 or higher
      * [MinGW](https://wiki.python.org/moin/WindowsCompilers#GCC_-_MinGW-w64_.28x86.2C_x64.29) 7 or higher
  * [Doxygen](https://www.doxygen.nl/download.html) 1.8.3 or higher
  * [Boost](https://www.boost.org/users/download)
  * [HDF5](https://www.hdfgroup.org/solutions/hdf5/) libraries
  
## Get the sources
For getting the sources files, just clone the github repository and create the build directory (out of the sources):

```sh
git clone https://github.com/gstlearn/gstlearn.git
cd gstlearn
```
Notes:
  * In the following, all instructions must be executed from a command prompt inside this *root* directory (thus the last command `cd gstlearn`)

## Library compilation & installation
For compiling and installing the *gstlearn* C++ Library, execute the following instructions (example given with static library (release version) and doxygen documentation generation):
#### GCC, Clang, MinGW, ...
```sh
cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Release
cmake --build build --target static
cmake --build build --target doxygen
sudo cmake --build build --target install
```
#### Microsoft Visual Studio, XCode, ...
```sh
cmake -Bbuild -H.
cmake --build build --target static --config Release
cmake --build build --target doxygen
sudo cmake --build build --target install # Not yet available !
```
Notes:
  * If you want to build and install the *Debug* version, you must replace `Release` by `Debug` above
  * The tests are compiled using the *shared* version of the *gstlearn* library
  * The *static* version of the library is mandatory for creating [pygstlearn package](https://github.com/gstlearn/pygstlearn)
  * You may need to precise the location of Boost or HDF5 header files, in that case, add the following in the first command above:
    * `-DBoost_INCLUDE_DIR="<path/to/boost/includes>"`
    * `-DHDF5_INCLUDE_DIR="<path/to/hdf5/includes>"`

## Usage
TODO: Instructions will come soon

## Required tools installation
### Linux (Ubuntu):
```sh
sudo apt install git
sudo apt install cmake
sudo apt install doxygen
sudo apt install libboost-all-dev
sudo apt install libhdf5-dev
```
Notes:
  * Under Linux, the GCC compiler and GNU make is already installed
  * If your Linux distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)

### MacOS:
```sh
brew install git
brew install cmake
brew install doxygen
brew install libboost-all-dev
brew install libhdf5-dev
```
Notes:
  * This is currently not tested - above packages may not exist
  * Under MacOS, the GCC (or Clang) compiler and GNU make is already installed
  * If your MacOS repository doesn't provide minimum required versions, please install manually (see provider website)
  
### Windows:
Download and install the following tools:
  * Git client [from here](https://gitforwindows.org) (Use default options during installation)
  * CMake tool [from here](https://cmake.org/download) (Check the 'Add CMake to the Path' option during installation)
  * Microsoft Visual C++ Compiler 14+ [from here](https://visualstudio.microsoft.com/visual-cpp-build-tools) (see Notes below) - OR - MinGW 7+ [from here](https://www.mingw-w64.org/downloads/)
  * Doxygen 1.8.3+ [from here](https://www.doxygen.nl/download.html) (Install in the directory *C:\\doxygen* for example)
  * Boost library [from here](https://www.boost.org/users/download) (Download and extract the zip file in *C:\\local\\* directory. If you choose another directory, CMake won't find it!)
  * HDF5 library [from here](https://www.hdfgroup.org/downloads/hdf5) (Download the pre-built binaries (zip), extract the zip file and execute the installer using default options)
    
Notes:
  * You must restart your computer after installing these requirements
  * The full Visual Studio C++ IDE is not necessary. You can 'only' download Visual Studio Build Tools (more details [here](https://stackoverflow.com/a/44398715)). Administrator rights are required.
  * If you prefer using another smaller compiler (i.e. MinGW), you could [try this](https://wiki.python.org/moin/WindowsCompilers#GCC_-_MinGW-w64_.28x86.2C_x64.29)
  * The *Path* environment variable must be updated to make *doxygen.exe* available (follow [this guide](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho) to add *C:\\doxygen\\bin* folder in the *Path* variable and restart Windows)

## Development
### Non-regression tests
### GCC, Clang, MinGW, ...
To launch non-regression tests, execute the following command:

```
cmake --build build --target build_test
cmake --build build --target test
```
#### Microsoft Visual Studio, XCode, ...
To build and launch non-regression tests, execute the following commands:

```
cmake --build build --target build_test
cd build
ctest -C Release
```
Notes:
  * If you want to run the *Debug* version of the tests, you must replace `Release` by `Debug` above (provided that *Debug* configuration has been built - see above)
  
### Clean
To clean (partially) the build, execute the following command:

```
cmake --build build --target clean
```
Note: If you really want to clean all files generated by CMake, you can remove *build* directory content by hand.

### Uninstall
To uninstall all the installed files (only the files, not the directories), execute this command:

```
sudo cmake --build build --target uninstall
```

### Generate the documentation
To (re)generate the documentation using doxygen, execute the command:

```
cmake --build build --target doxygen
```
The documentation is then available by opening the following HTML file with your favorite web-browser:

```
firefox build/doxygen/html/index.html
```

***

## License
MIT
2021 Team gstlearn

## Overview
*gstlearn* is the new cross-platform Geostatistics C++ Library proposed by MINES Paristech - PSL University. It offers to C++ users **all famous Geostatistical methodologies** developped and/or invented by the Geostatistic Team of the [Geosciences Research Center](https://www.geosciences.minesparis.psl.eu/)!<br/>
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
The *gstlearn* C++ Library is the direct successor of the Geoslib C/C++ Library which was proposed through the [RGeostats R package](http://cg.ensmp.fr/rgeostats)

## Requirements
This package has been successfully tested with Ubuntu 16.04 LTS, Ubuntu 18.04 LTS and Windows 10
For compiling and installing *gstlearn* C++ Library, the following tools must be available (See required tools installation instructions below):
  * [Doxygen](https://www.doxygen.nl/download.html) 1.8.3 or higher
  * [GCC](https://gcc.gnu.org) compiler 5.4 or higher (Linux/MacOS) or [Microsoft Visual C++ Compiler](https://visualstudio.microsoft.com/visual-cpp-build-tools) 14 or higher (Windows)
  * [Git](https://git-scm.com/downloads) client
  * [Boost](https://www.boost.org/users/download) library
  * [CMake](https://cmake.org/download) tool
  
## Get the sources
For getting the sources files, just clone the github repository:

```sh
git clone https://github.com/gstlearn/gstlearn.git
cd gstlearn
```

Notes:
  * In the following, all instructions must be executed from the gstlearn sources directory.
  
## Library compilation & installation
For compiling and installing the *gstlearn* library, execute the following instructions in a command prompt.

### Compile and install for all users
This requires Administrator rights:

```sh
cd build ; cmake  -DCMAKE_BUILD_TYPE=Release .. ; cd ..
cmake --build build --config Release
sudo cmake --build build --target install
```

### Compile and install for you only
You can choose a different installation directory by modifying the first CMake command above:
  
```
cd build ; cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=/home/fors/Programmes/gstlearn .. ; cd ..
cmake --build build --config Release
cmake --build build --target install
```

Additionnal options:
You may want to modify `make` behavior. In that case, add `--` at the end of the `cmake --build` command and complete:
  * If you want to use N CPU for compiling, add `-j N` after `--`
  * If you want to activate verbose mode, add `--no-print-directory VERBOSE=1` after `--`
  
Notes:
  * If you want to build and install the *Debug* version, you must replace *Release* by *Debug* above.
TODO : Ensure debug/release instructions are cross-platform

## Usage
TODO: Instructions will come soon

## Required tools installation

### Linux (Ubuntu):
Under Linux, the GCC compiler is already installed

```sh
sudo apt install doxygen
sudo apt install git
sudo apt install cmake
sudo apt install libboost-dev
```

### MacOS:
Under MacOS, the GCC (or Clang) compiler is already installed (Not tested)

```sh
brew install doxygen
brew install git
brew install cmake
brew install libboost-dev
```

### Windows:
Download and install the following tools:
  * Doxygen 1.8.3+ [from here](https://www.doxygen.nl/download.html) (installed in the directory *C:\\doxygen* for example)
  * Microsoft Visual C++ Compiler 14+ [from here](https://visualstudio.microsoft.com/visual-cpp-build-tools) (see Notes below)
  * Git client [from here](https://gitforwindows.org) (Use default options during installation)
  * Boost library [from here](https://www.boost.org/users/download) (download and extract the zip file anywhere)
  * CMake tool [from here](https://cmake.org/download)
    
Notes:
  * The full Visual Studio C++ IDE is not necessary. You can 'only' download Visual Studio Build Tools (more details [here](https://stackoverflow.com/a/44398715)). Administrator rights are required. If you prefer using another smaller compiler (i.e. MinGW), you could [try this](https://wiki.python.org/moin/WindowsCompilers#GCC_-_MinGW-w64_.28x86.2C_x64.29) (not tested)
  * The *Path* environment variable must be updated to make *doxygen.exe* available in the batch command line (follow [this guide](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho) to add *C:\\doxygen\\bin* folder in the *Path* variable and restart Windows)

## Development

### Clean & uninstall
To clean (partially) the build, execute the following command:

```
cmake --build build --target clean
```

To clean all CMake output, execute the following command:

```
rm -rf build/*
```

To uninstall all the installed files (only the files), execute this command:

```
cmake --build build --target uninstall
```



***

## License
MIT
2021 Team gstlearn

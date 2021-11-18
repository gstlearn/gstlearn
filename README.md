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
The *gstlearn* C++ Library is the direct successor of the Geoslib C/C++ Library which was proposed through the [RGeostats R package](http://cg.ensmp.fr/rgeostats)

## Requirements
This package has been successfully tested with Ubuntu 16.04 LTS, Ubuntu 18.04 LTS and Windows 10
For compiling and installing *gstlearn* C++ Library, the following tools must be available (See required tools installation instructions below):
  * [Git](https://git-scm.com/downloads) client
  * [CMake](https://cmake.org/download) tool 3.15 or higher
  * [GCC](https://gcc.gnu.org) compiler 5.4 or higher (Linux/MacOS) or [Microsoft Visual C++ Compiler](https://visualstudio.microsoft.com/visual-cpp-build-tools) 14 or higher (Windows) or MinGW 7 or higher (Windows also)
  * [Doxygen](https://www.doxygen.nl/download.html) 1.8.3 or higher
  * [Boost](https://www.boost.org/users/download) library
  
## Get the sources
For getting the sources files, just clone the github repository:

```sh
git clone https://github.com/gstlearn/gstlearn.git
cd gstlearn
```

Notes:
  * In the following, all instructions must be executed from the gstlearn sources directory
  
## Library compilation & installation
For compiling and installing the *gstlearn* library, execute the following instructions in a command prompt.

### Microsoft Visual Studio
(or any other multi-configuration generators)

```sh
cmake -Bbuild -H.
cmake --build build --config Release
sudo cmake --build build --target install
```

### Other compilers
(or any other single-configuration generators)

```sh
cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Release
cmake --build build
sudo cmake --build build --target install
```

Notes:
  * If you want to build and install the *Debug* version, you must replace `Release` by `Debug` above
  * If you don't want to generate doxygen documentation, use `-DEXCLUDE_DOXYGEN=1` in the first command
  * You may need to precise the location of Boost header files, in that case, add `-DBoost_INCLUDE_DIR="<path/to/boost_incs>"` in the first command

Additional options:
You may want to modify `make` behavior when running `cmake --build` command:
  * If you want to use N CPU for compiling, add `-j N` at the end
  * If you want to activate verbose mode, add `--no-print-directory VERBOSE=1` (Linux) or `--verbose` (Windows) at the end

## Usage
TODO: Instructions will come soon

## Required tools installation

### Linux (Ubuntu):
Under Linux, the GCC compiler is already installed.

```sh
sudo apt install git
sudo apt install cmake
sudo apt install doxygen
sudo apt install libboost-dev
```

Notes:
  * If your Linux distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)

### MacOS:
Under MacOS, the GCC (or Clang) compiler is already installed (Not tested)

```sh
brew install git
brew install cmake
brew install doxygen
brew install libboost-dev
```

Notes:
  * If your MacOS repository doesn't provide minimum required versions, please install manually (see provider website)
  
### Windows:
Download and install the following tools:
  * Git client [from here](https://gitforwindows.org) (Use default options during installation)
  * CMake tool [from here](https://cmake.org/download) (Check the 'Add CMake to the Path' option during installation)
  * Microsoft Visual C++ Compiler 14+ [from here](https://visualstudio.microsoft.com/visual-cpp-build-tools) (see Notes below) - OR - MinGW 7+ [from here] (https://www.mingw-w64.org/downloads/)
  * Doxygen 1.8.3+ [from here](https://www.doxygen.nl/download.html) (Install in the directory *C:\\doxygen* for example)
  * Boost library [from here](https://www.boost.org/users/download) (Download and extract the zip file in *C:\\local\\* directory. If you choose another directory, CMake won't find it!)
    
Notes:
  * The full Visual Studio C++ IDE is not necessary. You can 'only' download Visual Studio Build Tools (more details [here](https://stackoverflow.com/a/44398715)). Administrator rights are required.
  * If you prefer using another smaller compiler (i.e. MinGW), you could [try this](https://wiki.python.org/moin/WindowsCompilers#GCC_-_MinGW-w64_.28x86.2C_x64.29)
  * The *Path* environment variable must be updated to make *doxygen.exe* available in the batch command line (follow [this guide](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho) to add *C:\\doxygen\\bin* folder in the *Path* variable and restart Windows)

## Development
### Clean & uninstall
To clean (partially) the build, execute the following command:

```
cmake --build build --target clean
```
If you want to clean all CMake output, you can remove build directory:

```
rm -rf build
```

### Non-regression tests
#### Microsoft Visual Studio
To launch non-regression tests, execute the following command (from he *build* directory):

```
cd build
ctest -C Release
```
Notes:
  * If you want to run the *Debug* version of the tests, you must replace `Release` by `Debug` above (provided that *Debug* configuration has been built - see above)
  
#### Other compilers
To launch non-regression tests, execute the following command:

```
sudo cmake --build build --target test
```

### Uninstall
To uninstall all the installed files (only the files, not the directories), execute this command:

```
sudo cmake --build build --target uninstall
```

### Generate the documentation
To generate the documentation using doxygen, execute the command:

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

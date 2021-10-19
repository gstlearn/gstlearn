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
For using (compiling) *gstlearn* C++ Library, the following tools must be available (See required tools installation instructions below):
  * [Doxygen](https://www.doxygen.nl/download.html) 1.8.3 or higher
  * [GCC](https://gcc.gnu.org) compiler 5.4 or higher (Linux/MacOS) or [Microsoft Visual C++ Compiler](https://visualstudio.microsoft.com/visual-cpp-build-tools) 14 or higher (Windows)
  * [Git](https://git-scm.com/downloads) client
  * [Boost](https://www.boost.org/users/download) library (MacOS and Windows)
  * Following environment variables must be set:
      + ARCH to one of this: {linux64, macosx, windows} and
      + WIN_TYPE (Windows only) to one of {32, 64}
      + BOOST_DIR (MacOS and Windows only) to the Boost library root folder (See below)
  
## Library compilation
Cloning the repository and compiling (Currently, **only the Linux version** has been tested)
```sh
git clone https://github.com/gstlearn/gstlearn.git
cd gstlearn
make gstlearn
```

## Usage
TODO: Instructions will come soon

## Required tools installation
This package has been successfully tested with Ubuntu 18.04 LTS and Windows 10

### Linux (Ubuntu):
Under Linux, the GCC compiler is already installed
```sh
sudo apt install doxygen
sudo apt install git
```

### MacOS:
Under MacOS, the GCC (or Clang) compiler is already installed (Not tested)
```sh
brew install doxygen
brew install git
```

### Windows:
Download and install the following tools:
  * Doxygen 1.8.3+ [from here](https://www.doxygen.nl/download.html) (installed in the directory *C:\\doxygen* for example)
  * Microsoft Visual C++ Compiler 14+ [from here](https://visualstudio.microsoft.com/visual-cpp-build-tools) (see Notes below)
  * Git client [from here](https://gitforwindows.org)
  * Boost Library [from here](https://www.boost.org/users/download) (download and extract the zip file anywhere)
  
Notes:
  * The full Visual Studio C++ IDE is not necessary. You can 'only' download Visual Studio Build Tools (more details [here](https://stackoverflow.com/a/44398715)). Administrator rights are required. If you prefer using another smaller compiler (i.e. MinGW), you could [try this](https://wiki.python.org/moin/WindowsCompilers#GCC_-_MinGW-w64_.28x86.2C_x64.29) (not tested)
  * The *Path* environment variable must be updated to make *doxygen.exe* available in the batch command line (follow [this guide](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho) to add *C:\\doxygen\\bin* folder in the *Path* variable and restart Windows)

## Development
TODO: Instructions will come soon

***

## License
MIT
2021 Team gstlearn

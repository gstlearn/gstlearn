## Overview

The R *gstlearn* package is a cross-platform R Package wrapping the [gstlearn C++ Library](https://github.com/gstlearn/gstlearn). It offers to R users **all famous Geostatistical methodologies** developped and/or invented by the Geostatistic Team of the [Geosciences Research Center](https://www.geosciences.minesparis.psl.eu/)!<br/>
Copyright (c) MINES Paristech / PSL University

## References

The R *gstlearn* package is a R wrapper of the [gstlearn C++ Library](https://github.com/gstlearn/gstlearn).

## Requirements

For using this package, the requirements for building *gstlearn C++ library* must be installed:

* See [instructions here](https://github.com/gstlearn/gstlearn)
  
The following tools must be also available (See [required tools installation](#required-tools-installation) instructions below):

* [R](https://cran.r-project.org) 4 or higher with ggplot2 and gridExtra package
* [SWIG](http://www.swig.org/download.html) 4 or higher
* [RTools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html) for Windows user
  
## Build

### GCC, Clang, MinGW, ...

    cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Release
    cmake --build build --target r_build

or for those who prefer a single command line

    make r_build
  
### Microsoft Visual Studio, XCode, ...

    cmake -Bbuild -H.
    cmake --build build --target r_build --config Release

  
### Important Notes

Notes:

* If you want to build and install the *Debug* version, you must replace `Release` by `Debug` above
* You may need to precise the location of Boost or HDF5 installation directory (which contain *include* and *lib* folders). In that case, add the following in the first command above:
  * `-DBoost_ROOT=<path/to/boost>`
  * `-DHDF5_ROOT=<path/to/hdf5>`

## Usage

    # Load gstlearn package
    suppressWarnings(suppressMessages(library(gstlearn)))    
    # Grid size
    nx = 60
    ny = 30
    mygrid = DbGrid_create(c(nx,ny))
    # Add a uniform random field
    var = ut_vector_simulate_uniform(nx * ny)
    mygrid$addColumns(var, "var1")
    # Display the field
    plot.grid(mygrid, "var1")

## Changelog

Please, look at [CHANGES file](https://github.com/gstlearn/gstlearn/blob/main/CHANGES).

## Required tools installation

This package has been successfully tested with Ubuntu 16/18/20 LTS and Windows 10 (MacOS: not tested)

### Linux (Ubuntu):

Execute the following commands:

    sudo apt install r-base

Notes:

* If your Linux distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)

### MacOS:

Execute the following commands (Not tested):

    brew install r

Notes:

* If your MacOS distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)

### Windows:

Download and install the following tools:

* R [from here](https://cran.r-project.org)
* SWIG 4+ [from here](http://www.swig.org/download.html) (extract the archive in a directory of yours, let's say *C:\\swigwin-4.0.2*, see Notes below)
  
Notes:

* The *Path* environment variable must be updated to make *swig.exe* available in the batch command line (follow [this guide](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho) to add *C:\\swigwin-4.0.2* folder in the *Path* variable and restart Windows)
* The Windows C++ Compiler must be [RTools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html)

## Remove installed package

TODO


## Documentation

TODO

---

## License

MIT
2022 Team gstlearn

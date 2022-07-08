## Overview

The Python *gstlearn* package is a cross-platform Python Package wrapping the [gstlearn C++ Library](https://github.com/gstlearn/gstlearn). It offers to Python users **all famous Geostatistical methodologies** developped and/or invented by the Geostatistic Team of the [Geosciences Research Center](https://www.geosciences.minesparis.psl.eu/)!<br/>
Copyright (c) MINES Paris / PSL University

If you need to plot *gstlearn* outputs, you can import *gstlearn.plot* module which is based on *matplotlib* (see [modules](https://github.com/gstlearn/gstlearn/tree/main/python/modules) directory).

Some tutorials (Jupyter Notebooks) are provided in the [demo](https://github.com/gstlearn/gstlearn/tree/main/python/demo) directory.

Some tests (Python Scripts) are available in the [tests](https://github.com/gstlearn/gstlearn/tree/main/python/tests) directory.
  
## References

The Python *gstlearn* package is a Python wrapper of the [gstlearn C++ Library](https://github.com/gstlearn/gstlearn).

This package contains a copy of [doxy2swig](https://github.com/m7thon/doxy2swig) python script (see LICENSE.doxy2swig in *doc* folder).

## Requirements

1. For using this package, the requirements for building *gstlearn C++ library* must be installed: See [instructions here](https://github.com/gstlearn/gstlearn#required-tools-installation)
  
2. The following tools must be also available (See [required tools installation](#required-tools-installation) instructions below):

    * SWIG 4 or higher
    * Python 3 or higher with *pip*, *numpy*, *pybind11* and *matplotlib* modules installed
    * Optionnaly, following Python modules can also be installed [optional]: *pypandoc*, *geopandas*, *jupyter*

3. Finally, the source code of [gstlearn repository must be cloned](https://github.com/gstlearn/gstlearn#get-the-sources)

Note:

* In all commands below, users must use the correct `python` command (python command prompt) (can be `python3`, `python` or something else depending on your configuration)
  
## Installation
  
  These instructions will compile and install the Python package in your usual Python site-packages directory:
  
### Microsoft Visual Studio, XCode, ...

    cmake -Bbuild -H.
    cmake --build build --target python_doc --config Release
    cmake --build build --target python_install --config Release

### Important Notes

* Under Windows, using Visual Studio is mandatory for compiling Python packages
* The Python package documentation (python_doc target) is optional and requires Doxygen (see (here)[https://github.com/gstlearn/gstlearn#required-tools-installation])
* Using Visual Studio on a Windows where MingGW is also installed may need to add `-G "Visual Studio 16 2019"` in the first command (adapt version).
* If you want to build and install the *Debug* version, you must replace `Release` by `Debug` above
* You may need to precise the location of Boost or HDF5 installation directory (which contain *include* and *lib* folders). In that case, add the following in the first command above:
  * `-DBoost_ROOT=<path/to/boost>`
  * `-DHDF5_ROOT=<path/to/hdf5>`

## Usage

Simply import the *gstlearn* package and its plot module, then enjoy:

    # Import packages
    import numpy as np
    import gstlearn as gl
    import gstlearn.plot as gp
    # Grid size
    nx = 60
    ny = 30
    mygrid = gl.DbGrid.create([nx,ny],[1,1])
    # Add a gaussian random field
    var = np.random.randn(nx * ny)
    mygrid.addColumns(var, "var1")
    # Display the field
    gp.grid(mygrid, title="Gaussian random field", end_plot = True)

## Changelog

Please, look at [CHANGES file](https://github.com/gstlearn/gstlearn/blob/main/CHANGES).

## Required tools installation

This package has been successfully tested with Ubuntu 16/18/20 LTS and Windows 10 (MacOS: not tested)

### Linux (Ubuntu):

Execute the following commands:

    sudo apt install python3
    sudo apt install python3-pip
    sudo apt install swig
    python3 -m ensurepip --upgrade
    python3 -m pip install pybind11 numpy matplotlib
    python3 -m pip install pypandoc geopandas jupyter

Notes:

* If your Linux distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)
* According your Linux distribution you may have to replace pybind11 by the quoted string "pybind11[global]"

### MacOS:

Execute the following commands (Not tested):

    brew install python3
    brew install swig
    python3 -m ensurepip --upgrade
    python3 -m pip install pybind11 numpy matplotlib
    python3 -m pip install pypandoc geopandas jupyter

Notes:

* If your MacOS distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)
* According your MacOS distribution you may have to replace pybind11 by the quoted string "pybind11[global]"

### Windows:

Download and install the following tools:

* Python 3+ [from here](https://www.python.org/downloads) (which comes with pip and which is installed in *C:\\Python39* for example)
* SWIG 4+ [from here](http://www.swig.org/download.html) (extract the archive in a directory of yours, let's say *C:\\swigwin-4.0.2*, see Notes below)
* Pybind11, numpy and matplolib python modules by running following instructions in a command prompt:

    python -m pip install "pybind11[global]" numpy matplotlib
    python -m pip install pypandoc geopandas jupyter
  
Notes:

* The *Path* environment variable must be updated to make *swig.exe* (and *python.exe*) available in the batch command line (follow [this guide](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho) to add *C:\\swigwin-4.0.2* and *C:\\Python39* folder in the *Path* variable and restart Windows)
* The Windows C++ Compiler used must be the same that the one used for compiling Python (Visual C++). Using another compiler than Visual C++ is not supported.

## Remove installed package

To uninstall the Python Package, execute following command:

    python3 -m pip uninstall gstlearn

Note : You may need to directly modify your site-packages folder by:

* Removing the reference to the old gstlearn package version (see [this topic](https://stackoverflow.com/questions/43177200/assertionerror-egg-link-does-not-match-installed-location-of-reviewboard-at))
* Removing a line which contains gstlearn in the ./easy-install.pth file of the site-packages folder
* Removing all directories starting with *'~stlearn* from the site-packages folder

## Documentation

The classes and functions documentation is provided with the gstlearn package as html files generated by Doxygen [optional]. Please refer to gstlearn [README file](https://github.com/gstlearn/gstlearn) for more details. This documentation is also available using the `help` python command (if generated).

---

## License

MIT
2022 Team gstlearn

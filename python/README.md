## Overview
The *pygstlearn* package is a cross-platform Python Package wrapping the [gstlearn C++ Library](https://github.com/gstlearn/gstlearn). It offers to Python users **all famous Geostatistical methodologies** developped and/or invented by the Geostatistic Team of the [Geosciences Research Center](https://www.geosciences.minesparis.psl.eu/)!<br/>
Copyright (c) MINES Paristech / PSL University

If you need to plot *gstlearn* outputs, you can import *gstlearn.plot* module which is based on *matplotlib*.

Some tutorials (jupyter notebooks) are provided in the *demo* directory.

Some tests (Python scripts) are available in the *tests* directory.
  
## References
The *pygstlearn* package is the wrapper of the [gstlearn C++ Library](https://github.com/gstlearn/gstlearn).

This package contains a copy of [doxy2swig](https://github.com/m7thon/doxy2swig) python script (see LICENSE.doxy2swig in *swig* folder).

## Requirements
For using this package, the *gstlearn C++ library* must be installed:
  * See [instructions here](https://github.com/gstlearn/gstlearn)
  
The following tools must be also available (See [required tools installation](#required-tools-installation) instructions below):
  * [Python](https://www.python.org/downloads) 3 or higher with *pip*, *pybind11*, *matplotlib* and *scikit-build* modules installed
  * [SWIG](http://www.swig.org/download.html) 4 or higher

Note:
  * In all commands below, users must use the correct `python` command (python command prompt) (can be `python3`, `python` or something else depending on your configuration)
  
## Package installation from PyPi
TODO: Not yet available

## Package installation from sources
Cloning both repositories in the same folder and compiling/installing the Python package (this last command takes few minutes). Be sure not to forget the dot (at the end of the `pip install` command)
```sh
git clone https://github.com/gstlearn/pygstlearn.git
cd pygstlearn
python3 -m pip install .
```
Note:
  * The last command takes a while. You can add `-v` at the end to activate verbose mode.

## Usage
Simply import the package and its plot module, then enjoy:
```python
# Import packages
import numpy as np
import gstlearn as gl
import gstlearn.plot as gp
# Grid size
nx = 60
ny = 30
mygrid = gl.Db.createFromGrid([nx,ny],[1,1])
# Add a gaussian random field
var = np.random.randn(nx * ny)
mygrid.addFields(var, "var1")
# Display the field
gp.grid(mygrid, title="Gaussian random field", end_plot = True)
```

## Changelog
Please, look at [CHANGES file](https://github.com/gstlearn/gstlearn/blob/main/CHANGES) under gstlearn repository.

## Required tools installation
This package has been successfully tested with Ubuntu 18.04 LTS and Windows 10

### Linux (Ubuntu):
Execute the following commands:
```sh
sudo apt install python3
sudo apt install python3-pip
sudo apt install swig
python3 -m ensurepip --upgrade
python3 -m pip install pybind11 numpy matplotlib scikit-build
```
Notes:
* If your Linux distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)

### MacOS:
Execute the following commands (Not tested):
```sh
brew install python3
brew install swig
python3 -m ensurepip --upgrade
python3 -m pip install pybind11 numpy matplotlib scikit-build
```
Notes:
* If your MacOS distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)

### Windows:
Download and install the following tools:
  * Python 3+ [from here](https://www.python.org/downloads) (which comes with pip and which is installed in *C:\\Python39* for example)
  * SWIG 4+ [from here](http://www.swig.org/download.html) (extract the archive in a directory of yours, let's say *C:\\swigwin-4.0.2*, see Notes below)
  * Pybind11, numpy, matplolib and scikit-build python modules by running following instructions in a command prompt:
```
python -m pip install pybind11[global] numpy matplotlib scikit-build
```
  
Notes:
  * The *Path* environment variable must be updated to make *swig.exe* (and *python.exe*) available in the batch command line (follow [this guide](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho) to add *C:\\swigwin-4.0.2* and *C:\\Python39* folder in the *Path* variable and restart Windows)
  * The Windows C++ Compiler used by must be the same that the one used for compiling Python (Visual C++). If you prefer using another smaller compiler (i.e. MinGW), you could [try this](https://wiki.python.org/moin/WindowsCompilers#GCC_-_MinGW-w64_.28x86.2C_x64.29)

## Remove installed package
To uninstall the package, execute following command:
```
python3 -m pip uninstall gstlearn
```
Note : You may need to directly modify your site-packages folder by:
* Removing the reference to the old gstlearn package version (see [this topic](https://stackoverflow.com/questions/43177200/assertionerror-egg-link-does-not-match-installed-location-of-reviewboard-at))
* Removing a line which contains gstlearn in the ./easy-install.pth file of the site-packages folder
* Removing all directories starting with *'~stlearn* from the site-packages folder

## Documentation
The classes and functions documentation is provided with the gstlearn package as html files generated by doxygen. Please refer to gstlearn [README file](https://github.com/gstlearn/gstlearn) for more details. This documentation is also available using the `help` python command in the case that the documentation has been built (see above).
***

## License
MIT
2021 Team gstlearn

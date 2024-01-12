## Overview

**gstlearn** is the new cross-platform Geostatistics C++ library proposed by [MINES Paris - PSL University](https://www.minesparis.psl.eu/). It offers to users **all famous Geostatistical methodologies** developed and/or invented by the Geostatistic Team of the [Geosciences Research Center](https://www.geosciences.minesparis.psl.eu/).<br/>

The name 'gstlearn' stands for several purposes:

* GeoSTatistics & Machine Learning Library
* Geostatistical Spatio-Temporal Learning
* Learning Geosciences & Spatio-Temporal Models

*gstlearn* comes in various forms:

* A C++ library
* A Python package
* A R package

If you only want to use Python or R packages, you should switch to corresponding README:
* Python: [README](https://github.com/gstlearn/gstlearn/tree/main/python) or directly [install the package](https://soft.mines-paristech.fr/gstlearn/courses-latest/python/01_gstlearn_start.html)
* R: [README](https://github.com/gstlearn/gstlearn/tree/main/r) or directly [install the package](https://soft.mines-paristech.fr/gstlearn/courses-latest/r/01_gstlearn_start.html)

See [https://gstlearn.org](https://gstlearn.org) for more details.

## References

The *gstlearn* C++ library is the successor of the *Geoslib* C/C++ library which was proposed through the [RGeostats R package](http://cg.ensmp.fr/rgeostats).

The *gstlearn* C++ library is developed by the [Geostatistics Group](https://www.geosciences.minesparis.psl.eu/en/presentation/geostatistics) of the Geosciences Center ([MINES Paris - PSL University](https://www.minesparis.psl.eu) - France)

The *gstlearn* C++ library :
* is a derivative work based on the [*swigex0* project](https://github.com/fabien-ors/swigex0).
* depends on third-party libraries and source codes (see [below](#credits)).
* comes with several data files that are used for our [documentation](https://gstlearn.org/?page_id=50) (tutorials and courses).

See [credits](#credits) below.

## How to cite

When using the *gstlearn* C++ Library, please, use the following to cite us in any publication or results for which **gstlearn** has been used:

```
--------------------------------------------------
gstlearn
Geostatistics and Machine Learning toolbox
Copyright © MINES Paris - PSL University
Free download from https://gstlearn.org
--------------------------------------------------
```

You may be interested in the citation file [gstlearn.bib](https://soft.mines-paristech.fr/gstlearn/gstlearn.bib)

## Report a bug

To report a bug or contact us:
* feel free to post an [issue](https://github.com/gstlearn/gstlearn/issues) or
* visit our [Help](https://gstlearn.org/?page_id=468) page

## Development

### Requirements

This library has been successfully tested with Ubuntu 18/20/22 LTS and Windows 10 (MacOS: not tested).
For **compiling and installing** *gstlearn* C++ library, the following tools must be available (See [required tools installation](#required-tools-installation) instructions below):

* Git client 2.30 or higher
* CMake tool 3.20 or higher
* A C++ compiler among:
  * Linux/MacOS:
    * GCC 5.4 or higher
    * Clang 12 or higher (not tested)
  * Windows:
    * Python users: Microsoft Visual Studio C++ 14 or higher
    * R users: MinGW 7 (RTools 4.2) or higher
* Boost header files 1.65 or higher
* Doxygen [Optional] 1.8.3 or higher with LaTeX and Ghostscripts
* HDF5 [Optional] C++ library and header files 1.8 or higher

See [required tools installation](#required-tools-installation) instructions below

### Get the sources

For getting the sources files, just clone the github repository:

```
git clone https://github.com/gstlearn/gstlearn.git
cd gstlearn
```

Next time, you will only need to pull the repository (If you have some local undesirable modifications, you have to revert them and execute the pull, otherwise do not execute `git reset`):

````
cd gstlearn
git reset --hard
git pull
````

### C++ Library Compilation & Installation

For compiling and installing the *gstlearn* C++ shared library, execute the following instructions from the *root* directory of *gstlearn*. Please note that you can choose another destination folder (currently named *build*).

#### GCC, Clang, MinGW, ...

...or any other single configuration compiler:

```
cmake -Bbuild -S. -DCMAKE_BUILD_TYPE=Release
cmake --build build --target shared
cmake --build build --target install
```

or even faster:

```
make
```

#### Microsoft Visual Studio, ...

...or any other multiple configurations compiler:

```
cmake -Bbuild -S.
cmake --build build --target shared --config Release
cmake --build build --target install --config Release
```

### Execute Non-regression Tests

The `check*` targets bring some required runtime customization, so do not use the standard *ctest* command for triggering the non-regression tests.

To build and launch non-regression tests, execute the following commands:

#### GCC, Clang, MinGW, ...

...or any other single configuration compiler:

```
cmake --build build --target build_tests
cmake --build build --target check_cpp
cmake --build build --target check_data
```

or even faster:

```
make check_cpp
make check_data
```

#### Microsoft Visual Studio, ...

...or any other multiple configurations compiler:

```
cmake --build build --target build_tests --config Release
cmake --build build --target check_cpp --config Release
cmake --build build --target check_data --config Release
```

### Clean

To clean (partially) the build, execute the following command:

```
cmake --build build --target clean
```

Notes:

* If you really want to clean all files generated by CMake, you can remove *build* directory content by hand. Linux, MacOS or MinGW users may use the `clean_all` target from the shortcut Makefile inside the top level directory:

````
make clean_all
````

### Usage

Please, look at *tests* [C++ code](https://github.com/gstlearn/gstlearn/tree/main/tests) in order to learn how to use the *gstlearn* C++ library. You can generate the source code documentation using [Doxygen](#generate-the-documentation).

### Required Tools Installation

These tools are needed for compiling the *gstlearn* C++ library. Please note that HDF5 and Doxygen (and Latex) installation are optional.

If you modified your system, you must reinstall the requirements from scratch following next instructions. You must delete 'gstlearn' existing source folders (if so).

Note :

* In case of issues, see [Important Notes below](#important-notes).

#### Linux (Ubuntu)

```
sudo apt install git
sudo apt install cmake
sudo apt install texlive-latex-recommended
sudo apt install texlive-science
sudo apt install doxygen
sudo apt install libboost-all-dev
sudo apt install libhdf5-dev
```

#### MacOS

```
brew install git
brew install cmake
brew install texlive-latex-recommended
brew install texlive-science
brew install doxygen
brew install libboost-all-dev
brew install libhdf5-dev
```

Notes:

* These instructions for MacOS are currently not tested - above packages may not exist
  
#### Windows - Microsoft Visual Studio

These requirements are also recommended to people who wants to compile *gstlearn* Python package. If you want to compile *gstlearn* R package under Windows, you should look at the next section.

##### Install all tools

Download and install the following tools using default options during installation:

1. Git client [from here](https://gitforwindows.org) (*Setup program* [exe])
2. CMake tool [from here](https://cmake.org/download) (*Windows Installer* [msi], check the *'Add CMake to the system PATH for all users'* option during installation)
3. Microsoft Visual Studio (Community) [from here](https://visualstudio.microsoft.com/fr/thank-you-downloading-visual-studio/?sku=Community&channel=Release&version=VS2022&source=VSLandingPage&cid=2030&passive=false) (*VisualStudioSetup.exe* - only select the *Visual Studio Desktop C++* component)
4. Boost library [from here](https://www.boost.org/users/download) (*Archive file* [zip] to be extracted in a folder of your choice, but not in the *gstlearn* source folder - and remind that folder)
5. HDF5 library (optional) [from here](https://www.hdfgroup.org/downloads/hdf5) (*Pre-built binaries* [zip] to be extracted, then, execute the *installer* [msi] - and remind the installation folder)
6. Doxygen (optional) [from here](https://www.doxygen.nl/download.html) (*Binary distribution* [setup.exe] - remind the installation folder, we assume it is `C:\Program Files\doxygen`)
7. LaTeX and Ghostscripts following instructions [here](https://www.doxygen.nl/manual/install.html#install_bin_windows)

##### Update the Path environment variable

The *Path* environment variable (*System variables*) must be updated to make *doxygen.exe* available in the batch command line:

1. Follow [this guide](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10) to add *bin* directory from the *Doxygen* installation folder in the *Path* System variable (i.e: `C:\Program Files\doxygen\bin`)
2. Restart Windows


#### Windows - MinGW (via RTools):

These requirements are also recommended to people who wants to compile *gstlearn* R package. If you want to compile *gstlearn* Python package under Windows, you should look at the previous section. This is not the only way to install MinGW. But using MinGW provided with RTools permits us to also handle *gstlearn* R package compilation.

##### Install R and RTools

Remove all R and RTools previous installation and download and install the following tools using default options:

1. R [from here](https://cran.r-project.org/bin/windows/base) (*Setup program* [exe] - remind the installation folder, assume it is `C:\Program Files\R\R-4.2.2`)
2. RTools [from here](https://cran.r-project.org/bin/windows/Rtools) (RTools *Installer* [exe] - remind the installation folder, assume it is `C:\rtools42`)

Notes:

* Choose the corresponding RTools version according to the R version installed
* Instructions in this section are **valid since R v4.2** (for older versions please contact us)
* RTools is not the unique way to install MinGW on Windows, but it is our preferred way as we can handle R packages compilation

##### Update the Path environment variable

The *Path* environment variable (*System variables*) must be updated to make *R.exe* available in the batch command line:

1. Follow [this guide](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10) to add `bin` directory from the *R* installation folder in the *Path* variable (ie: `C:\Program Files\R\R-4.2.2\bin`)
2. Restart Windows

##### Add MSYS2 Required Packages

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
pacman -Sy mingw-w64-x86_64-texlive-latex-recommended
pacman -Sy mingw-w64-x86_64-texlive-science
pacman -Sy mingw-w64-x86_64-doxygen
````

### Important Notes

* If your system distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)
* If you plan to generate the documentation, add `-DBUILD_DOXYGEN=ON` to the first cmake command above.
* If you don't know how to execute github commands or you experience a 'password authentication' problem, you may [read this](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token).
* Currently, **HDF5 is not supported** when compiling *gstlearn* C++ library **under Windows and MacOS**. *gstlearn* won't link against HDF5 and GibbsMMulti::setFlagStoreInternal(false) feature won't be available.
* The default installation directory named *gstlearn_install* is located in your *Home*. If you want to change it, you can add `-DCMAKE_INSTALL_PREFIX="path/of/gstlearn/install/dir"` to the first cmake command above.
* If you don't want HDF5 support, add `-DUSE_HDF5=OFF` to the first cmake command above. If you use the shortcut Makefile, you can use `USE_HDF5=0` after the `make` command
* Only the *shared* library (built by default) is installed when compiling *gstlearn* C++ library. If you want to compile the *static* version, you must replace *shared* by *static* target above.
* Using MinGW on a Windows where another compiler is also installed may need to add `-G "MSYS Makefiles"` in the first cmake command above.
* Using Visual Studio on a Windows where another compiler is also installed may need to add `-G "Visual Studio 16 2019"` in the first command (adapt version).
* If you want to build and install the *Debug* version, you must replace `Release` by `Debug` above. If you use the shortcut Makefile, you can use `DEBUG=1` after the `make` command
* You may need to precise the location of Boost, HDF5 or Doxygen installation directory. In that case, add the following variables in the first cmake command above:
   * `-DBoost_ROOT="path/to/boost"`
   * `-DHDF5_ROOT="path/to/hdf5"`
   * `-DDoxygen_ROOT="path/to/doxygen"`

### Uninstall the Library

To uninstall all the installed files (only the files, not the directories), execute this command:

```
cmake --build build --target uninstall
```

or faster:

```
make uninstall
```

### Generate the Documentation

The Doxygen HTML documentation is optional (not included in the installation by default). If you want to generate it, execute the command:

```
cmake -Bbuild -S. -DBUILD_DOXYGEN=ON
cmake --build build --target doxygen
```

or faster (for Makefile user):

```
make doxygen
```

The documentation is then available by opening the following HTML file with your favorite web-browser:

```
firefox build/doxygen/html/index.html
```

## Credits

### Derivative work

The *gstlearn* C++ library is a derivative work based on the *swigex0* project (see license in [**doc/licenses**](doc/licenses)):

| Name           | License        | URL                                                            | Copyright
|----------------|----------------|----------------------------------------------------------------|-----------
| swigex0        | MIT            | https://github.com/fabien-ors/swigex0                          | Copyright 2023, Fabien Ors


### Third party libraries

The *gstlearn* C++ library depends on the following third-party source code, slightly modified and compiled in separate libraries (see [**3rd-party**](3rd-party) folder):

| Name           | License        | URL                                                            | Copyright
|----------------|----------------|----------------------------------------------------------------|-----------
| csparse        | LGPL v2.1      | https://people.math.sc.edu/Burkardt/c_src/csparse/csparse.html | Copyright 2006, Timothy A. Davis
| stripack (GMT) | LGPL v3        | https://www.generic-mapping-tools.org                          | Copyright(c) 2020, the GMT Team


### Third party source code

The *gstlearn* C++ library includes external source codes (see licenses notices in [**doc/licenses**](doc/licenses)):

| Name                    | License        | URL                                            | Copyright
|-------------------------|----------------|------------------------------------------------|-----------
| clustering              | Python License | http://bonsai.hgc.jp/~mdehoon/software/cluster | Copyright (C) 2002 Michiel Jan Laurens de Hoon
| fft                     | see licenses   | https://netlib.org/go/fft-olesen.tar.gz        | Copyright(c)1995,97 Mark Olesen
| ball                    | MIT            | https://42.fr                                  | Copyright(c) 2017 Eung Bum Lee
| sparseinv (SuiteSparse) | BSD 3-clause   | http://www.suitesparse.com                     | Copyright 2011, Timothy A. Davis
| vtk (VisIt)             | BSD 3-clause   | https://visit.llnl.gov                         | Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC

The *gstlearn* C++ library also depends on the following third-party libraries (see licenses notices in [**doc/licenses**](doc/licenses)):

| Name           | License        | URL                                                            | Copyright
|----------------|----------------|----------------------------------------------------------------|-----------
| Boost          | see licenses   | https://www.boost.org                                          | see Boost headers
| HDF5           | see licenses   | https://www.hdfgroup.org                                       | Copyright 2006 by The HDF Group


### Data files

The *gstlearn* C++ library comes with several data files that are used for our [documentation](https://gstlearn.org/?page_id=50) (tutorials and courses).

Here are the credits and licenses for the different data files available in each directories of [**doc/data**](doc/data):

| Name          | License        | URL                                            | Copyright
|---------------|----------------|------------------------------------------------|-----------
| Alluvial      | Etalab v 2.0   | https://infoterre.brgm.fr                      | Copyright (C) Sept. 2021, BRGM
| benchmark     | see (1)        | http://www.ai-geostats.org                     | Copyright (C) 2003, Dubois, G., Malczewski, J. and De Cort, M.
| boundaries    | N/A            | https://www.naturalearthdata.com               | Made with Natural Earth
| BRGM          | see (2)        | https://brgm.fr                                | Copyright (C) Sept. 2021, BRGM derived from Bardossy et al. 1999
| Chamaya       | CC BY 4.0      | https://gstlearn.org                           | Copyright (C) 2010, Renard, D. & Beucher, H.
| FKA           | CC BY 4.0      | https://gstlearn.org                           | Copyright (C) 1999, Renard, D.
| halieutic     | CC BY 4.0      | https://ices-library.figshare.com              | Copyright (C) 2017, ICES see (3)
| Meshings      | CC0 1.0 Univ.  | https://www.cs.cmu.edu/~kmcrane                | Copyright (C) 2020, Crane, Keenan and Pinkall, Ulrich and Schröder, Peter
| PluriGaussian | CC BY 4.0      | https://gstlearn.org                           | Copyright (C) 1999, Renard, D.
| Pollution     | CC BY 4.0      | http://rgeostats.free.fr                       | Copyright (C) 2000, Team RGeostats
| Scotland      | CC BY 4.0      | http://rgeostats.free.fr                       | Copyright (C) 2000, Team RGeostats
| Selectivity   | CC BY 4.0      | https://gstlearn.org                           | Copyright (C) 2010, Freulon, X.

(1) Reference: Dubois, G., Malczewski, J. and De Cort, M. (2003). Mapping radioactivity in the environment. Spatial Interpolation Comparison 1997 (Eds.). EUR 20667 EN, EC. 268 p.

(2) Bárdossy, A., H. Giese, J. Grimm-Strele, and K. Barufke (2003), SIMIK + GIS—Implementierte Interpolation von Grundwasserparametern mit Hilfe von Landnutzungs- und Geologiedaten. Hydrol. Wasserwirt., 47, 13–20.

(3) ICES. 2017. Handbook of Geostatistics in R for Fisheries and Marine Ecology. ICES COOPERATIVE RESEARCH REPORT. Vol. 338, 184 pp. https://doi.org/10.17895/ices.pub.3717

---

## License

*gstlearn* C++ Library is distributed under the license:

BSD 3-clause

2024 Team gstlearn

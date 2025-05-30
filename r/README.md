## gstlearn: The Geostatistics &amp; Machine Learning R Package
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13343742.svg)](https://doi.org/10.5281/zenodo.13343742)

The **gstlearn** R package is a cross-platform R package wrapping the [gstlearn C++ Library](https://gstlearn.org). It offers to R users **all famous Geostatistical methodologies** developed and/or invented by the Geostatistic Team of the [Geosciences Research Center](https://www.geosciences.minesparis.psl.eu/)!

All *gstlearn* outputs can be ploted using *plot* generic function which relies on *ggplot* package.


## References

The *gstlearn* R package is a R wrapper of the [gstlearn C++ Library](https://gstlearn.org). It's the successor of the [RGeostats R package](http://cg.ensmp.fr/rgeostats).

This package contains a modified copy of [findR.cmake](https://github.com/root-project/root) script (see LICENSE.root in *root* folder).q

The *gstlearn* R package is a derivative work based on the *swigex0* project: [https://github.com/fabien-ors/swigex0](https://github.com/fabien-ors/swigex0)


## How to cite

When using the *gstlearn* R Package, please, use this to cite us in any publication or results for which **gstlearn** has been used:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13343742.svg)](https://doi.org/10.5281/zenodo.13343742)

You may be interested in the citation file [gstlearn.bib](https://soft.mines-paristech.fr/gstlearn/gstlearn.bib)

## Installation

For using this R Package you only need R 4.2 (or higher) and execute the following R command:

```
install.packages("gstlearn",repos="https://soft.mines-paristech.fr/cran")
```

Note: With slow Internet connection, you may need to increase the default timeout (60s) (gstlearn package is around 17Mo!)

```
options(timeout=1000)
```

## Usage

We refer the reader to this [course page](https://soft.mines-paristech.fr/gstlearn/courses-latest/r/01_gstlearn_start.html) for an introduction and important information about R gstlearn package.

Simply load the *gstlearn* and ggplot2 R package, then enjoy:

```
# Load gstlearn package
library(gstlearn)
# Grid size
nx = 60
ny = 30
mygrid = DbGrid_create(c(nx,ny), c(1,1))
# Add a gaussian random field
var = rnorm(nx * ny)
mygrid$addColumns(var, "var1", ELoc_Z())
# Display the field
plot.init(asp=1) + plot.raster(mygrid) + plot.decoration(title="Gaussian random field")
```

Some tutorials (R Markdown) are provided in the *demo* directory [here](https://github.com/gstlearn/gstlearn/tree/main/doc/demo/r) and their HTML rendering is provided [here](https://soft.mines-paristech.fr/gstlearn/demos-latest/r/).

Some tests (R Scripts) are available in the [tests](https://github.com/gstlearn/gstlearn/tree/main/tests/r) directory of the *gstlearn* github repository.


## Changelog

Please, look at [CHANGES file](https://github.com/gstlearn/gstlearn/blob/main/CHANGES).


## Developments

### Requirements

For building the *gstlearn* R package, the requirements for building *gstlearn C++ library* must be installed beforehand. Then, the following additional tools must be also available:

* SWIG 4.2.0 **customized by Fabien Ors** (not the official version!)
* R 4.2 or higher
* RTools 4.2 or higher (for Windows users only)
* *roxygen2* R package (for generating R documentation)
* *ggplot2*, *vctrs*, *ggpubr*, *ggrepel*, *ggnewscale* R packages [Optional] (only for plotting)
* *FNN*, *rgl*, *lares*, *Matrix*, *knitr*, *tidyr*, *geigen* and *callr* R packages [Optional] (only for testing R Markdown scripts)

If you modified your system (or if you installed a new version or RTools), you must reinstall the requirements from scratch following next instructions. You must delete 'gstlearn' and 'swig' existing source folders (if so).

Note :

* In case of issues, see [Important Notes below](#important-notes).

#### Linux (Ubuntu)

1. Install *gstlearn* C++ library requirements for Linux [here](https://github.com/gstlearn/gstlearn#linux-ubuntu)

2. Remove any previous installation of SWIG (if any)

3. Then, execute the following commands:

````
sudo apt install r-base
sudo apt install bison
sudo apt install pcre2-devel # Ubuntu 18
sudo apt install libpcre2-dev # Ubuntu 20
````

4. In a directory of your choice, get the source of SWIG 4.2.0 [customized] by executing following commands (in the same shell):

````
git clone https://github.com/fabien-ors/swig.git
cd swig
````

Next time (if a new version of SWIG 4.2.0 [customized] is delivered), you will only need to pull the repository

````
cd swig
git pull
````

5. Then compile and install SWIG 4.2.0 [customized] :

````
cmake -Bbuild
cd build
make
sudo make install
````

6. Finally, install the R packages from an R command prompt (only roxygen2 is mandatory):

````
install.packages(c("roxygen2"), repos="https://cloud.r-project.org")
install.packages(c("ggplot2", "vctrs", "ggpubr", "ggrepel", "ggnewscale"), repos="https://cloud.r-project.org")
install.packages(c("FNN", "rgl", "lares", "Matrix", "knitr", "callr", "tidyr", "geigen"), repos="https://cloud.r-project.org")
````

#### MacOS

1. Install *gstlearn C++ library* requirements for MacOS [here](https://github.com/gstlearn/gstlearn#macos)

2. Remove any previous installation of SWIG (if any)

3. Then, execute the following commands (Not tested):

````
brew install r
brew install bison
brew install pcre2-devel
````

4. In a directory of your choice, get the source of SWIG 4.2.0 [customized] by executing following commands (in the same shell):

````
git clone https://github.com/fabien-ors/swig.git
cd swig
````

Next time (if a new version of SWIG 4.2.0 [customized] is delivered), you will only need to pull the repository

````
cd swig
git pull
````

5. Then compile and install SWIG 4.2.0 [customized] :

````
cmake -Bbuild
cd build
make
sudo make install
````

6. Finally, install the R optional packages from an R command prompt (only roxygen2 is mandatory):

````
install.packages(c("roxygen2"), repos="https://cloud.r-project.org")
install.packages(c("ggplot2", "vctrs", "ggpubr", "ggrepel", "ggnewscale"), repos="https://cloud.r-project.org")
install.packages(c("FNN", "rgl", "lares", "Matrix", "knitr", "callr", "tidyr", "geigen"), repos="https://cloud.r-project.org")
````

Note :

* These instructions for MacOS are currently not tested - above packages may not exist

#### Windows

1. Install *gstlearn C++ library* requirements for Windows (RTools) [here](https://github.com/gstlearn/gstlearn#windows---mingw-via-rtools)

2. Remove any previous installation of SWIG (if any)

3. Launch *mingw64.exe* in RTools installation directory (i.e.: `C:\rtools42`) and pin the icon to the task bar

````
pacman -Sy bison
pacman -Sy mingw-w64-x86_64-pcre2
````

4. In a directory of your choice, get the source of SWIG 4.2.0 [customized] by executing following commands (in the same shell):

````
git clone https://github.com/fabien-ors/swig.git
cd swig
````

Next time (if a new version of SWIG 4.2.0 [customized] is delivered), you will only need to pull the repository

````
cd swig
git pull
````

5. Then compile and install SWIG 4.2.0 [customized] :

````
cmake -G "MSYS Makefiles" -Bbuild -DCMAKE_INSTALL_PREFIX:PATH=/mingw64/
cd build
make
make install
````

6. Finally, install the R optional packages from an R command prompt (only roxygen2 is mandatory):

````
install.packages(c("roxygen2"), repos="https://cloud.r-project.org")
install.packages(c("ggplot2", "vctrs", "ggpubr", "ggrepel", "ggnewscale"), repos="https://cloud.r-project.org")
install.packages(c("FNN", "rgl", "lares", "Matrix", "knitr", "callr", "tidyr", "geigen"), repos="https://cloud.r-project.org")
````

### Installation from Source

1. For getting the *gstlearn* R package sources files, just clone the github repository:

````
git clone https://github.com/gstlearn/gstlearn.git
cd gstlearn
````

Next time, you will only need to pull the repository (If you have some local undesirable modifications, you have to revert them and execute the pull, otherwise do not execute `git reset`):

````
cd gstlearn
git reset --hard
git pull
````

2. Then, these instructions will compile and install the *gstlearn* R package in your usual R libraries directory:

````
cmake -Bbuild -S. -DBUILD_R=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target r_install
````

or even faster:

```
make r_install
```

### Execute Non-regression Tests

The `check*` targets bring some required runtime customization, so do not use the standard *ctest* command for triggering the non-regression tests.

To build and launch non-regression R tests, you need to execute the following command:

```
cmake --build build --target check_r
```

or even faster:

```
make check_r
```

### Important Notes

* ggnewscale and ggplot2 packages must be installed/updated together! Otherwise, you could have troubles when using plotting feature.
* Under Linux or MacOS, if you don't have sudo permissions, you may have to install swig in a folder of your choice. In that case, use `-DCMAKE_INSTALL_PREFIX:PATH=/home/user/Programs/swig4.2.0b` (adapt installation folder) in the `cmake` command above.
* If your system distribution repository doesn't provide minimum required versions, please install the tools manually (see provider website)
* You may need to reconnect to your session after installing some requirements
* If you experience the following issue: `Error: ERROR: no permission to install to directory...`, we suggest you to run the `install.packages` command (at least one time). This will create a *personal R library folder* having writing permissions.
* If you plan to generate the documentation, add `-DBUILD_DOC=ON` to the first cmake command above. Then users will be able to execute `make doxygen`.
* If you don't know how to execute github commands or you experience a 'password authentication' problem, you may [read this](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token).
* Under Windows, using RTools is mandatory for compiling R packages
* Under Windows, you may need to add `-G "MSYS Makefiles"` to the first cmake command above
* If you want to build and install the *Debug* version, you must replace `Release` by `Debug` above
* You may need to precise the location of Boost, Eigen3, SWIG, Doxygen, HDF5 or NLopt installation directory. In that case, add the following variables in the first cmake command above:
  * `-DBoost_ROOT="path/to/boost"`
  * `-DEigen3_ROOT="path/to/eigen3"`
  * `-DSWIG_ROOT="path/to/swig"`
  * `-DDoxygen_ROOT="path/to/doxygen"`
  * `-DHDF5_ROOT="path/to/hdf5"`
  * `-DNLopt_ROOT="path/to/nlopt"`

### Remove Installed Package

Execute the following command from an R command prompt:

```
remove.packages("gstlearn")
```

---

## License

GPL v3

2025 Team gstlearn

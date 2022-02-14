## Windows compilation of gstlearn using Rtools40

### Install & customize Rtools40

#### Download and install

* Download rtools40 for Windows 64 bits (rtools40-x86_64.exe) from here : https://cran.r-project.org/bin/windows/Rtools/rtools40.html
* Install rtools40 by executing the downloaded program
* Add the following directories to *Path* environment variable (see [this guide](https://stackoverflow.com/questions/44272416/how-to-add-a-folder-to-path-environment-variable-in-windows-10-with-screensho):
  + %RTOOLS40_HOME%\\usr\\bin
  + %RTOOLS40_HOME%\\mingw64\\bin
* Restart Windows

#### Add/upgrade rtools
* From a Windows command prompt, execute following instructions:

```sh
#pacman -Syu
#pacman -S mingw-w64-x86_64-toolchain
pacman -S mingw-w64-x86_64-cmake
pacman -S mingw-w64-x86_64-hdf5
pacman -S mingw-w64-x86_64-boost
```

### Then compile HDF5 (obsolete)

The idea is to:
- use MSYS generator from RTools,
- build C++ static library and prevent from building fortran/java libraries and tests

#### Following these steps 

1- Download HDF5 CMake version here : https://www.hdfgroup.org/downloads/hdf5/source-code/#

2- Extract and enter directory CMake-hdf5-1.12.1 (adapt version number)

3- Remove all CRLF in characters strings from src/H5EInit.h file !!! Example :

Line 26, replace :

if((msg = H5E__create_msg(cls, H5E_MAJOR, "Invalid arguments to routine
"))==NULL)

by 

if((msg = H5E__create_msg(cls, H5E_MAJOR, "Invalid arguments to routine"))==NULL)

4- Then execute the following (adapt version number and installation directory below):

mkdir build
cd build
cmake -G "MSYS Makefiles" -DHDF5_GENERATE_HEADERS:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=OFF -DDEFAULT_API_VERSION:STRING=v110 -DCMAKE_BUILD_TYPE:STRING=Release -DHDF5_BUILD_FORTRAN:BOOL=OFF -DHDF5_BUILD_CPP_LIB:BOOL=ON -DHDF5_BUILD_JAVA:BOOL=OFF -DCMAKE_INSTALL_PREFIX:PATH=C:\local -DHDF5_ENABLE_SZIP_SUPPORT:BOOL=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DHDF5_BUILD_TOOLS:BOOL=OFF ..\hdf5-1.12.1
make
make install

Notes:
  - Currently (Jan. 2022), the gstlearn code is based on HDF5 1.10. (So we must use the flag -DDEFAULT_API_VERSION:STRING=v110)
  - We need static c++ library only (thus -DBUILD_SHARED_LIBS:BOOL=OFF & -DHDF5_BUILD_CPP_LIB:BOOL=ON)
  - We don't regenerate headers because that fails (-DHDF5_GENERATE_HEADERS:BOOL=OFF)
  - We don't need tools, tests, fortran, java and szip/z_lib support (I bet !)  


### Finally, compile gstlearn & RGeostats (obsolete)

cd gstearn
cmake -Bbuild_msys -H. -G "MSYS Makefiles" -DHDF5_ROOT:PATH=C:\local -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_C_FLAGS="-DH5_USE_110_API"
cmake --build build_msys --target static

cd ../rgeostats
R CMD INSTALL --preclean --no-multiarch --build --latex --example --html RGeostats

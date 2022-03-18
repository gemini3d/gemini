# Build Gemini3D with Intel oneAPI on MacOS

Intel oneAPI does not come with MPI on MacOS.
The gemini3d/external repo will build MPI for oneAPI on MacOS.

## setup Intel oneAPI

oneAPI relies on the underlying AppleClang compiler for libc and libstdc++.
Xcode is required:

```sh
xcode-select --install
```

Install
[oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)
with these options:

* Math Kernel Library (oneMKL)

Install
[oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html)
with these options:

* Intel C++ compiler
* Intel Fortran compiler

We recommend making a little shell script named like "~/intel_oneapi.sh".
The contents of this script would be like:

```sh
source /opt/intel/oneapi/setvars.sh

export LIBRARY_PATH=$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
export CPLUS_INCLUDE_PATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include

export CXXFLAGS=-I$CPLUS_INCLUDE_PATH
export CFLAGS=$CXXFLAGS

export FC=ifort CC=icc CXX=icpc
```

To enable oneAPI in this Terminal:

```sh
source ~/intel_oneapi.sh
```

## One-time Gemini3D external library setup

```sh
git clone https://github.com/gemini3d/external

cmake -S external -B external/build -DCMAKE_INSTALL_PREFIX=~/lib_intel

cmake --build external/build
```

that installs Gemini3d external libraries (including MPI) under ~/lib_intel.
This path is arbitrary but should be distinct between compilers.

## Build and Test Gemini3D

```sh
git clone https://github.com/gemini3d/gemini3d

cmake -S gemini3d -B build/gemini3d -G Ninja -DCMAKE_PREFIX_PATH=~/lib_intel

cmake --build gemini3d/build

ctest --test-dir gemini3d/build
```
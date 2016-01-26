# 2D Time Domain BEM

## Introduction
The 2D Time Domain Boundary Element Method (TDBEM) method is an electromagnetic simulation scheme capable of accurate and dispersion free description of outwardly radiating fields.

This implementation in C++ computes the 2D TDBEM operators, and is intended to be used in collaboration with my Matlab program, found in my [BEUT] repository.

The project repository has the following directory structure:
* 2DTDBEM - contains the license, readme and install script
 * 2DTDBEM - contains the source code and CMakeLists file
 * build - where the built binaries go
 * CMake - additional tools for CMake
 * input - where the input Matlab files go
 * results - where the output Matlab files go

## Installation
To install, use the CMakeLists file. There are dependencies on the following libraries:
* OpenMP (optional but recommended for parallel computing and faster run-times)
* Zlib (optional but recommended for compression)
* HDF5 (optional but recommended for large files)
* MatIO (required)
* Armadillo (required)

### Linux
From a fresh install (e.g. Amazon Web Services EC2 Ubuntu), download the project to a custom directory, then run the `install.sh` script:
```
sudo apt-get install git
git clone https://github.com/dan-phd/2DTDBEM.git
cd 2DTDBEM
chmod +x install.sh
sudo ./install.sh
```
Alternatively, you can install the above libraries (in order) yourself using the standard `./configure` and `sudo make install` commands.

If you don't have root privileges, you will need to install the libraries in a local folder, and use `./configure --prefix=/home/<local_lib_folder>` or `cmake . -DCMAKE_INSTALL_PREFIX=/home/<local_lib_folder>`.

If HDF5 is installed, MatIO should be configured using:
```
./configure --with-default-file-ver=7.3
```
Once the libraries have been installed, 2DTDBEM is installed using:
```
cd 2DTDBEM/build
sudo cmake ..
sudo make install
cd ..
```
The program options can then be viewed with `./bin/2DTDBEM`.

If you get `error wile loading shared libraries`, use `export LD_LIBRARY_PATH=/usr/local/lib/` for root users, or `export LD_LIBRARY_PATH=/home/<local_lib_folder>/lib` otherwise.

### Windows
First download and unzip the armadillo and MatIO libraries to a low level directory (such as `C:\build`).

Edit the user environment variables for your machine to include the directories to the unzipped libraries in variables named `ARMADILLO_ROOT` and `MATIO_ROOT`.

To run the code using Visual Studio in Windows, open the `2DTDBEM.sln` file located in the top directory.

Make sure that the correct solution platform is being used (Win32 or x64) and check the following settings are configured correctly:
* Configuration Properties -\> C/C++ -\> General -\> Additional Include Directories:
 * `$(ARMADILLO_ROOT)\include`
 * `$(MATIO_ROOT)\include`
* Configuration Properties -\> Linker -\> General -\> Additional Library Directories:
 * `$(ARMADILLO_ROOT)\include\examples\lib_win64`
 * `$(MATIO_ROOT)\visual_studio\x64\Debug` (or `Release` depending on solution configuration)
* Configuration Properties -\> Linker -\> Input -\> Additional Dependencies:
 * `lapack_win64_MT.lib`
 * `blas_win64_MT.lib`
 * `libmatio.lib`

After building the Visual Studio project but before running, make sure that `lapack_win64_MT.dll`, `blas_win64_MT.dll`, and `libmatio.dll` is copied to the runtime directory i.e. the same directory as the `2DTDBEM.exe file`.

To run the application straight from Visual Studio, you can append program arguments in Configuration Properties -\> Command Arguments. The program
can then be run without the debugger by pressing `Ctrl` + `F5`.

## Beginning test
Once you are all set up, run the following initial test to check everything works:
```
./bin/2DTDBEM --test computeConvolutions
```


[BEUT]: https://github.com/dan-phd/BEUT
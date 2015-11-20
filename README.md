# 2D Time Domain BEM

## Introduction
The 2D Time Domain Boundary Element Method (TDBEM) method is an electromagnetic simulation scheme capable of accurate and dispersion free description of outwardly radiating fields.

This implementation in C++ computes the 2D TDBEM operators, and is intended to be used in collaboration with my Matlab program, found in my [BEUT] repository.

## Installation
To install, use the CMakeLists file. There are dependencies on the following libraries:
* OpenMP (optional but recommended for faster run-times)
* Zlib
* HDF5 (optional but recommended for large files)
* MatIO
* Armadillo

### Linux
From a fresh install (e.g. Amazon Web Services EC2 Ubuntu), download the project to a custom directory, then run the `install.sh` script:
```
sudo apt-get install git
git clone https://github.com/dan-phd/2DTDBEM.git
cd 2DTDBEM
chmod u+x install.sh
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


## Beginning test
Once you are all set up, run the following initial test to check everything works:
```
./bin/2DTDBEM --test computeConvolutions
```


[BEUT]: https://github.com/dan-phd/BEUT
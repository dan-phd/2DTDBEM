# 2D Time Domain BEM

## Introduction
The 2D Time Domain Boundary Element Method (TDBEM) method is an electromagnetic simulation scheme capable of accurate and dispersion free description of outwardly radiating fields.

This implementation in C++ computes the 2D TDBEM operators, and is intended to be used in collaboration with my Matlab program, found in my [BEUT] repository.

## Installation
To install, use the CMakeLists file. There are dependencies on the following libraries:
* Armadillo
* MatIO
* Zlib
* HDF5 (optional but recommended for large files)
* OpenMP (optional but recommended for faster run-times)

### Linux
From a fresh install (e.g. Amazon Web Services EC2 Ubuntu):
```
sudo apt-get update
sudo apt-get install build-essential liblapack-dev libarpack++2-dev libopenblas-dev cmake libhdf5-dev zlib1g-dev -y
mkdir build && cd build
wget ftp://ftp.hdfgroup.org/HDF5/current/src/hdf5-*.tar.gz
wget http://sourceforge.net/projects/matio/files/latest/download
wget http://sourceforge.net/projects/arma/files/armadillo-6.100.1.tar.gz
tar zxf hdf5….tar.gz
cd hdf5…
sudo ./configure
sudo make install
cd ..
tar zxf download
cd matio…
./configure --with-default-file-ver=7.3
sudo make install
cd ..
tar zxf armadillo….tar.gz
cd armadillo…
sudo ./configure
sudo make install
```
#### Installing TDBEM
Transfer/download all files in the `2DTDBEM` to a custom directory, move to this directory and type the following:
```
mkdir build && cd build
cmake ..
make install
cd ..
export LD_LIBRARY_PATH=/usr/local/lib/
mkdir results
./bin/2DTDBEM
```


## Beginning test
Once you are all set up, run the following initial test to check everything works:
```
./bin/2DTDBEM --test computeConvolutions
```


[BEUT]: https://github.com/dan-phd/BEUT
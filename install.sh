#!/usr/bin/env bash

ROOT_UID=0     # Only users with $UID 0 have root privileges.
LINES=50       # Default number of lines saved.
E_XCD=86       # Can't change directory?
E_NOTROOT=87   # Non-root exit error.


# Run as root, of course.
if [ "$UID" -ne "$ROOT_UID" ]
then
  echo "Must be root to run this script."
  exit $E_NOTROOT
fi

if [ -n "$1" ]
# Test whether command-line argument is present (non-empty).
then
  lines=$1
else  
  lines=$LINES # Default, if not specified on command-line.
fi

# get required packages
apt-get update
apt-get install build-essential liblapack-dev libarpack++2-dev libopenblas-dev cmake libhdf5-dev zlib1g-dev unzip -y

# download and install HDF5, MatIO and Armadillo
mkdir build && cd build
wget http://downloads.sourceforge.net/project/matio/matio/1.5.2/matio-1.5.2.zip http://sourceforge.net/projects/arma/files/armadillo-6.100.1.tar.gz
unzip matio-1.5.2.zip
cd matio-1.5.2 || {
  echo "Cannot change to necessary directory. Check dowload link." >&2
  exit $E_XCD;
}
./configure --with-default-file-ver=7.3
make install
cd ..
tar zxf armadillo-6.100.1.tar.gz
cd armadillo-6.100.1 || {
  echo "Cannot change to necessary directory. Check dowload link." >&2
  exit $E_XCD;
}
./configure
make install
cd ..

exit 0
#  A zero return value from the script upon exit indicates success
#+ to the shell.
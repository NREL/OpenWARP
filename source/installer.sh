#!/bin/bash
set -e
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in?page=1&tab=votes#tab-top
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

INSTALL_PATH="$DIR"

ROOT="$DIR"

echo "Installing to $ROOT"

echo "Installing GFortran and gcc"
sudo apt-get --yes --force-yes install build-essential gfortran gcc

echo "Installing cmake"
sudo apt-get --yes --force-yes install cmake

echo "Installing Blas and Lapack"
sudo apt-get --yes --force-yes install liblapack-dev libblas-dev


echo "Installing nglib"
sudo apt-get --yes --force-yes install libnglib-4.9.13 netgen netgen-headers libnglib-dev

echo "Installing OpenCASCADE Community Edition"
sudo apt-get --yes --force-yes install liboce-foundation-dev liboce-modeling-dev  liboce-ocaf-dev liboce-ocaf-lite-dev liboce-visualization-dev oce-draw

echo "Installing vtk"
sudo apt-get --yes --force-yes install libvtk5-dev

echo "Installing python dependencies"
sudo apt-get --yes --force-yes install python-numpy python-scipy python-matplotlib python-h5py python-cherrypy3 python-pip ipython ipython-notebook python-pandas python-sympy python-nose


mkdir -p /tmp/OpenWarp
rm -rf /tmp/OpenWarp/*

echo "Installing mesh generator nglib-mesh"

cd /tmp/OpenWarp/

cmake "$ROOT/openwarpgui/bundled/mesh-generator/src"
make
cp nglib_mesh "$ROOT/openwarpgui/bundled/mesh-generator/bin"


echo "Installing Nemoh Fortran"
rm -rf /tmp/OpenWarp/*

cmake -DCMAKE_Fortran_COMPILER="gfortran" "$ROOT/NemohImproved/Nemoh"
make

cp libnemoh.so "$ROOT/openwarpgui/bundled/simulation/libs"

export LD_LIBRARY_PATH="$ROOT/openwarpgui/bundled/simulation/libs"
export LDFLAGS="-L$ROOT/openwarpgui/bundled/simulation/libs"


cd "$ROOT/openwarpgui/nemoh"

chmod +x "$INSTALL_PATH"/openwarpgui/bundled/mesh-generator/bin/nglib_mesh

chmod +x "$INSTALL_PATH"/openwarpgui/bundled/paraview_linux/bin/paraview


# Avoiding any problem by using the system pip and python
sudo /usr/bin/pip install -r "$ROOT/openwarpgui/requirements.txt"
/usr/bin/python setup.py cleanall
/usr/bin/python setup.py build_ext --inplace

echo "OpenWarp Installation successfully completed"





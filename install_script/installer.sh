#!/bin/bash

OSTYPE=`uname`

#Stop shell scripts if commands are returing non-zero value
#set -e

# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in?page=1&tab=votes#tab-top


SOURCE="${BASH_SOURCE[0]}"

while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
SOURCE="$(readlink "$SOURCE")"
[[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done

DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
PARENTDIR="$(dirname "${DIR}")"
INSTALL_PATH="${PARENTDIR}/source"
ROOT="${PARENTDIR}/source"

NEMOH_FORTRAN="${ROOT}/NemohImproved/Nemoh"
FORTRAN_BUILD="${ROOT}/NemohImproved/FortranBuild"

echo "Root is" $ROOT

if [ "$OSTYPE" = "Linux" ];then
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

	echo "Installing HDF5-tools"
	sudo apt-get --yes --force-yes install hdf5-tools

	echo "Installing vtk"
	sudo apt-get --yes --force-yes install libvtk5-dev

	echo "Installing python dependencies"
	sudo apt-get --yes --force-yes install python-numpy python-scipy python-matplotlib python-h5py python-cherrypy3 python-pip ipython ipython-notebook python-pandas python-sympy python-nose python-progressbar python-vtk


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

elif [ "$OSTYPE"="Darwin" ];then
	
	USERNAME=`id -un`
	echo "Username is:"${USERNAME}
	
	echo "Installing command line tools for Xcode"
	xcode-select --install
	
	echo "Installing Homebrew"
	ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

	echo "Brew Version"
	brew --version	

	echo "Brew Doctor"
	brew doctor
	
	echo "Installing curl"
	brew install curl 

	echo "Installing GCC and GFortran"
	brew install gcc 	

	echo "Installing CMake"
	brew install cmake

	echo "Downloading Anaconda Command Line Installer"
	if ! [ -f Anaconda2-4.1.1-MacOSX-x86_64.sh ]; then
		curl -O http://repo.continuum.io/archive/Anaconda2-4.1.1-MacOSX-x86_64.sh
	else
		echo "Andaconda Command Line Installer is downloaded already"
	fi
	
	
	echo "Instaling Anaconda Command Line Installer"
	bash ./Anaconda2-4.1.1-MacOSX-x86_64.sh	
	
	echo "Adding Anaconda to the path"
	export PATH=/Users/${USERNAME}/anaconda2/bin:$PATH

	echo "Create new directory FORTRAN_BUILD"
	mkdir ${FORTRAN_BUILD}
	cd ${FORTRAN_BUILD}
	
	# It is important if we are creating installer
	echo "Make sure the dynamic version of quad-math library is not in the path"
	mv /usr/local/lib/gcc/4.9/libquadmath.0.dylib /usr/local/lib/gcc/4.9/disable_libquadmath.0.dylib
	mv /usr/local/lib/gcc/4.9/libquadmath.dylib /usr/local/lib/gcc/4.9/disable_libquadmath.dylib
	
	echo "Compiling Nemoh Fortran"
	cmake -DCMAKE_FortranE_COMPILER="gfortran" $NEMOH_FORTRAN
	make
	
	echo "copying libnemoh.dylib from FORTRAN_BUILD to lib/directory inside the Anaconda installation Root"
	cp $FORTRAN_BUILD/libnemoh.dylib /Users/${USERNAME}/anaconda2/lib
	
	echo "Downloading ParaView "
	cd ${DIR}
	# -L follows the redirect in cURL
	curl -L http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v4.1&type=binary&os=osx&downloadFile=ParaView-4.1.0-Darwin-64bit.dmg
	
	# Above failes than we can download from a dropbox location
	curl -L -o ParaView-4.1.0-Darwin-64bit.dmg  https://www.dropbox.com/s/ocdmgl6dwzaq60l/ParaView-4.1.0-Darwin-64bit.dmg?dl=1
	echo "Installing Parview"
	#bash ./ParaView-4.1.0-Darwin-64bit.dmg

	echo "Copy Paraview App Folder"
	# TODO:
	#Mount the DMG file
	hdiutil attach -mountpoint ${DIR}ParaView-4.1.0-Darwin-64bit.dmg
	# copy the Paraview.app file
	# copy from disk to src/openwarpgui/bundled
	sudo cp ParaView-4.1.0-Darwin-64bit.app ${ROOT}/openwarpgui/bundled/
	
	echo "Command to test the library path is correctly set "
	(test -e /Users/${USERNAME}/anaconda2/lib/libnemoh.dylib && echo ’Success’ ) || echo ’Error:Nemoh library not found’.
	
	echo "make sure the nglib-mesh directory is executable"
	chmod +x ${ROOT}/openwarp/bundled/meshgenerator/build/nglib-mesh
	
	echo "preventing dylib errors (if any)"
	sudo python ${ROOT}/openwarpgui/bundled/mesh-generator/lib/update_dylib.py

	echo "Installing pip"
	sudo easy_install pip

	echo "Installing the remaining python dependencies"
	sudo pip install -r ${ROOT}/openwarpgui/requirements.txt
	python ${ROOT}/openwarpgui/nemoh/setup.py cleanall
        python ${ROOT}/openwarpgui/nemoh/setup.py build_ext --inplace
	
	# To know the currently configured path 
	otool -L ${ROOT}/openwarpgui/nemoh/solver_fortran.so
	# install_name_tool -change $Current $New
	install_name_tool -change "libnemoh.dylib" "/Users/${USERNAME}/anaconda2/lib/libnemoh.dylib" ${ROOT}/openwarpgui/nemoh/solver_fortran.so
	
	echo "OpenWarp Installationg Succesfully Completed"	
	
else
 	echo "This OS type is not supported!  "
fi



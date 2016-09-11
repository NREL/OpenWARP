#!/bin/bash
set -e

OSTYPE=`uname`

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

cd "$INSTALL_PATH"/openwarpgui

chmod +x "$INSTALL_PATH"/openwarpgui/bundled/mesh-generator/bin/nglib_mesh

chmod +x "$INSTALL_PATH"/openwarpgui/bundled/paraview_linux/bin/paraview

export LD_LIBRARY_PATH="$ROOT/openwarpgui/bundled/simulation/libs"
export LDFLAGS="-L$ROOT/openwarpgui/bundled/simulation/libs"
export PYTOHNPATH="$INSTALL_PATH"/openwarpgui

mkdir -p ~/OpenWarpFiles/temp
mkdir -p ~/OpenWarpFiles/user_data

echo "Starting OpenWarp GUI Server"

# Avoiding any problem by using the system pip and python
# Run the main.py without sudo to avoid any problem with export
# If it is needed to run this as sudo, then the export commands should be run as sudo too
if [ "$OSTYPE"="Linux" ];then
	/usr/bin/python "$INSTALL_PATH/openwarpgui/main.py"
elif [ "$OSTYPE"="Darwin" ];then
	python "$INSTALL_PATH/openwarpgui/main.py"
fi

#!/usr/bin/env bash
# Source in bash before building/running this solver and its local ICSFoam library.
export DBNS_IHB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${FOAM_BASHRC:-/usr/lib/openfoam/openfoam2112/etc/bashrc}" || return $?
export ICSFOAM_BUILD_SRC="$DBNS_IHB_DIR/Make/ICSFoam"
export FOAM_USER_LIBBIN="$DBNS_IHB_DIR/platforms/$WM_OPTIONS/lib"
export FOAM_USER_APPBIN="$DBNS_IHB_DIR/platforms/$WM_OPTIONS/bin"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH"
export PATH="$FOAM_USER_APPBIN:$PATH"

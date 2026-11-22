#!/usr/bin/env bash
# Reuse the same v2112 ICSFoam library as the project's incompressible solver.
set -eo pipefail
: "${WM_PROJECT_DIR:?Source the OpenFOAM v2112 bashrc first}"
[[ $WM_PROJECT_VERSION == v2112 ]] || { echo "Requires OpenFOAM v2112" >&2; exit 1; }
here=$(cd "$(dirname "$0")" && pwd)
apps=$(dirname "$here")
inc=$apps/dbnsIncompressibleHBFoam
export ICSFOAM_BUILD_SRC=${ICSFOAM_BUILD_SRC:-$inc/Make/ICSFoam}
lib=${ICSFOAM_LIB_DIR:-$inc/platforms/$WM_OPTIONS/lib}
[[ -f $lib/libICSFoam.so && -f $ICSFOAM_BUILD_SRC/lnInclude/HBZone.H ]] || {
    echo "Build the project ICSFoam library first, or set ICSFOAM_LIB_DIR and ICSFOAM_BUILD_SRC." >&2
    exit 1
}
export FOAM_USER_LIBBIN=$lib
export FOAM_USER_APPBIN=$here/platforms/$WM_OPTIONS/bin
mkdir -p "$FOAM_USER_APPBIN"
build=$(mktemp -d)
echo "HB application build directory: $build"
mkdir "$build/HBacoustics"
cp "$here/runHBAcoustics.H" "$build/HBacoustics/"
for name in dbnsIncompressibleHBFoam dbnsFullyImplicitHBFoam; do
    app=$build/$name
    mkdir -p "$app/Make"
    cp "$apps/$name/"*.{C,H} "$app/"
    cp "$apps/$name/Make/"{files,options} "$app/Make/"
    sed -i 's|../../../src/lnInclude|$(ICSFOAM_BUILD_SRC)/lnInclude|' "$app/Make/options"
    wmake -j "${BUILD_JOBS:-2}" "$app"
done

#!/usr/bin/env bash
# Build a private v2112 libAcoustics with the HB wall-data adapter.
set -eo pipefail
: "${WM_PROJECT_DIR:?Source the OpenFOAM v2112 bashrc first}"
[[ $WM_PROJECT_VERSION == v2112 ]] || { echo "Requires OpenFOAM v2112" >&2; exit 1; }
here=$(cd "$(dirname "$0")" && pwd)
upstream=${1:?Pass the path to a local libAcoustics checkout}
dest=${2:?Pass a private library output directory}
mkdir -p "$dest"
export FOAM_USER_LIBBIN=$(cd "$dest" && pwd)
build=$(mktemp -d)
echo "libAcoustics build directory: $build"
revision=1ac0140e4dcfe8c8add29293c793fc57d73854ec
git -C "$upstream" archive "$revision" Sources/lib | tar -x -C "$build"
cd "$build"
find Sources/lib -type f \( -name '*.H' -o -name '*.C' -o -path '*/Make/*' \) -exec sed -i 's/\r$//' {} +
patch --fuzz=0 -p1 < "$here/libAcoustics-wall-interface.patch"
cp "$here/hbSolidWallFWH.C" Sources/lib/
sed -i '/^LIB =/i hbSolidWallFWH.C\n' Sources/lib/Make/files
if [[ ! -f ${FFTW_ARCH_PATH:-}/include/fftw3.h && -f /usr/include/fftw3.h ]]; then
    export FFTW_ARCH_PATH=/usr
fi
wmake -j "${BUILD_JOBS:-4}" libso Sources/lib

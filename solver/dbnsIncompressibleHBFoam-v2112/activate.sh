#!/usr/bin/env bash
# Load OpenFOAM-v2112 first, then source this file in Bash.
if [[ ! ${WM_PROJECT_DIR:-} || ! ${WM_OPTIONS:-} || ${WM_PROJECT_VERSION#v} != 2112 ]]; then
    echo "Load the OpenFOAM-v2112 environment before sourcing activate.sh." >&2
    return 1
fi
_hb_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P) || return 1
export PATH="$_hb_dir/platforms/$WM_OPTIONS/bin:$PATH"
export LD_LIBRARY_PATH="$_hb_dir/platforms/$WM_OPTIONS/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
unset _hb_dir

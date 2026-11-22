#!/bin/bash
# Python runs on the host; the stock v2112 image supplies OpenFOAM and g++.
set -euo pipefail
project_dir="$(cd "$(dirname "$0")/.." && pwd)"
cd "$project_dir"
python3 tests/generate_cases.py
docker run --rm --platform linux/amd64 --entrypoint /bin/bash \
    -v "$project_dir:/work" -w /work opencfd/openfoam-default:2112 -lc '
    source /usr/lib/openfoam/openfoam2112/etc/bashrc
    set -e
    export FOAM_USER_APPBIN=/work/bin
    mkdir -p "$FOAM_USER_APPBIN"
    ./Allwmake > tests/build.log 2>&1
    g++ -std=c++11 -O2 -Wall -Wextra -pedantic tests/algebra_test.cpp -o tests/algebra_test
    ./tests/algebra_test > tests/algebra.log
    for name in oscillatoryChannel amiCavity amiOscillatoryCavity; do
        bash "tests/cases/$name/Allmesh"
    done'
python3 tests/negative_cases.py
docker run --rm --platform linux/amd64 --entrypoint /bin/bash \
    -v "$project_dir:/work" -w /work opencfd/openfoam-default:2112 -lc '
    source /usr/lib/openfoam/openfoam2112/etc/bashrc
    set -e
    for name in oscillatoryChannel amiCavity amiOscillatoryCavity checkerboardChannel steadyAmiCavity; do
        /work/bin/fullyImplicitHBFoam -case "tests/cases/$name" > "tests/cases/$name/log.solver" 2>&1
    done
    set +e
    /work/bin/fullyImplicitHBFoam -case tests/cases/maxNewtonZero > tests/cases/maxNewtonZero/log.solver 2>&1
    status=$?
    set -e
    test "$status" -eq 2
    test ! -d tests/cases/maxNewtonZero/1'
for name in oscillatoryChannel amiCavity amiOscillatoryCavity checkerboardChannel steadyAmiCavity; do
    python3 tests/check_results.py "tests/cases/$name" --report "tests/cases/$name/verification.json"
done
python3 - <<'PY'
import json
from pathlib import Path
p = Path('tests/cases/maxNewtonZero')
assert 'FAILED: no converged solution written' in (p/'log.solver').read_text()
assert not (p/'1').exists()
(p/'verification.json').write_text(json.dumps({'passed': True, 'expectedExitCode': 2, 'observedExitCode': 2, 'solutionWritten': False}, indent=2)+'\n')
print('PASS: maxNewtonZero returns 2 and writes no solution')
PY

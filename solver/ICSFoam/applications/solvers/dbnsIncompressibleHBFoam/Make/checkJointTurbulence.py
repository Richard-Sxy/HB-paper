#!/usr/bin/env python3
"""One Newton step for uniform, shear-free k-epsilon (source env.sh first).

Physical residuals are -epsilon and -C2*epsilon**2/k. Their joint Newton
direction is (-k,-epsilon). The 20% relative-change guard therefore gives
k=0.8 and epsilon=0.08 from (1,0.1), even with HB, convection and MPI enabled.
Exit code 2 from the solver is expected: this is a step test, not a steady root.
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    assert binary.is_file(), binary
    source = Path(__file__).resolve().parent / "turbulenceCheck/perturbed"
    case = Path(tempfile.mkdtemp(prefix="ihb-joint-analytic-"))
    for name in ("constant", "system"):
        shutil.copytree(source / name, case / name)
    for rank in range(2):
        dest = case / f"processor{rank}"
        dest.mkdir()
        for name in ("constant", "0"):
            shutil.copytree(source / dest.name / name, dest / name)
    for path in case.rglob("fvOptions"):
        path.write_text("FoamFile { version 2.0; format ascii; "
                        "class dictionary; object fvOptions; }\n")
    for rank in range(2):
        for field, value in (("k", "1"), ("epsilon", "0.1"), ("nut", "0.9")):
            for path in (case / f"processor{rank}/0").rglob(field):
                path.write_text(re.sub(r"(uniform\s+)[-+0-9.eE]+",
                                       r"\g<1>" + value, path.read_text()))
    path = case / "system/controlDict"
    text = re.sub(r"\bendTime\s+[^;]+;", "endTime 0.1;", path.read_text())
    text = re.sub(r"\bwriteInterval\s+[^;]+;", "writeInterval 1;", text)
    path.write_text(text)
    path = case / "system/fvSolution"
    path.write_text(path.read_text().replace(
        "pseudoTime {", "pseudoTime { hbTurbulenceMethod phaseImplicit; "
        "phaseTurbulenceNewton true; phaseTurbulenceDirections 40; "
        "phaseTurbulenceCycles 5; phaseTurbulenceRelTol 1e-5;"))
    log = case / "log.check"
    with log.open("w") as stream:
        run = subprocess.run(["mpirun", "-np", "2", str(binary), "-parallel",
                              "-case", str(case)], stdout=stream,
                             stderr=subprocess.STDOUT, timeout=90)
    text = log.read_text()
    assert run.returncode == 2, (run.returncode, log)
    assert "IHBPhaseTurbulenceTransaction accepted 1" in text, log
    assert "IHBFinished converged=0" in text, log
    result = {"case": str(case), "stepTestOnly": True, "fields": {}}
    for field, expected in (("k", 0.8), ("epsilon", 0.08)):
        values = []
        paths = list(case.glob(f"processor*/0.1/subTimeLevel*/{field}"))
        assert len(paths) == 6, paths
        for path in paths:
            text = path.read_text()
            match = re.search(r"internalField\s+uniform\s+([^;]+);", text)
            if match:
                values.append(float(match[1]))
            else:
                match = re.search(r"internalField\s+nonuniform\s+List<scalar>"
                                  r"\s+\d+\s*\((.*?)\)", text, re.S)
                assert match, path
                values.extend(map(float, match[1].split()))
        error = max(abs(value - expected) for value in values)
        assert error < 1e-5, (field, error, case)
        result["fields"][field] = {"expected": expected, "maxError": error,
                                   "minimum": min(values), "maximum": max(values)}
    (case / "check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

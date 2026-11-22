#!/usr/bin/env python3
"""Compare ASCII OpenFOAM channel U against an independent discrete HB solution.

The reference is the finite-volume Stokes system solved in complex arithmetic,
not the solver's Newton iteration. Cell ordering is blockMesh x-fastest.
"""
import argparse
import json
from pathlib import Path
import re
import sys


def internal_vectors(path, n):
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)
    uniform = re.search(r"internalField\s+uniform\s*\(([^)]+)\)\s*;", text)
    if uniform:
        return [tuple(map(float, uniform.group(1).split()))] * n
    match = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"Cannot parse ASCII internalField in {path}")
    values = [tuple(map(float, s.split())) for s in re.findall(r"\(([^()]+)\)", match.group(2))]
    if len(values) != n or int(match.group(1)) != n or any(len(v) != 3 for v in values):
        raise ValueError(f"Unexpected vector count in {path}: expected {n}, got {len(values)}")
    return values


def check(case, time_name, tolerance):
    metadata = json.loads((case / "metadata.json").read_text())
    if metadata["case"] != "oscillatoryChannel":
        raise ValueError("The checker accepts only the oscillatoryChannel case")
    nx, ny = metadata["nx"], metadata["ny"]
    errors, transverse, continuum_errors = [], [], []
    for phase in metadata["phases"]:
        velocity = internal_vectors(case / time_name / phase["name"] / "U", nx*ny)
        error = max(abs(velocity[j*nx+i][0] - phase["discreteUbyY"][j]) for j in range(ny) for i in range(nx))
        cont_error = max(abs(velocity[j*nx+i][0] - phase["continuumUbyY"][j]) for j in range(ny) for i in range(nx))
        cross_error = max(max(abs(v[1]), abs(v[2])) for v in velocity)
        errors.append(error)
        transverse.append(cross_error)
        continuum_errors.append(cont_error)
        print(f"{phase['name']}: discrete max|Ux-error|={error:.8g}; max|Uy,Uz|={cross_error:.8g}; continuum error={cont_error:.8g}")
    result = {"passed": max(errors+transverse) <= tolerance,
              "time": time_name, "discreteMaxAbsError": max(errors), "transverseMaxAbsVelocity": max(transverse),
              "continuumMaxAbsError": max(continuum_errors), "acceptanceTolerance": tolerance,
              "note": "Continuum error includes spatial discretization and is reported, not used as the nonlinear-solve criterion."}
    print(json.dumps(result, indent=2))
    return result["passed"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--time", default="1", help="Written solution directory; default 1")
    parser.add_argument("--tolerance", type=float, default=2e-5)
    args = parser.parse_args()
    sys.exit(0 if check(args.case, args.time, args.tolerance) else 1)

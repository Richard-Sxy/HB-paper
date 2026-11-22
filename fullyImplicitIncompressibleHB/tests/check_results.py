#!/usr/bin/env python3
"""Check fresh solver logs plus written solutions for the generated HB cases."""
import argparse
import json
import math
from pathlib import Path
import re
import sys

from check_channel import internal_vectors


def internal_scalars(path, n):
    text = re.sub(r"/\*.*?\*/|//[^\n]*", "", path.read_text(), flags=re.S)
    uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if uniform:
        return [float(uniform.group(1))] * n
    match = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"Cannot parse scalar field {path}")
    values = list(map(float, match.group(2).split()))
    if len(values) != n or int(match.group(1)) != n:
        raise ValueError(f"Unexpected scalar count in {path}")
    return values


def verify(case, time_name="1"):
    metadata = json.loads((case / "metadata.json").read_text())
    log = (case / "log.solver").read_text()
    metrics = re.findall(r"maxMomentum\s+(\S+)\s+maxContinuity\s+(\S+)\s+maxAMIMass\s+(\S+)", log)
    if not metrics:
        raise ValueError("Missing final physical residual audit in log.solver")
    momentum, continuity, ami_mass = map(float, metrics[-1])
    momentum_fluxes = re.findall(r"maxAMIMomentum\s+(\S+)", log)
    ami_momentum = float(momentum_fluxes[-1]) if momentum_fluxes else math.inf
    gcl = [float(x) for x in re.findall(r"normalized GCL error\s+(\S+)", log)]
    coverage = [float(x) for x in re.findall(r"coverage/area error\s+(\S+)", log)]
    result = {"case": metadata["case"], "time": time_name,
              "maxMomentum": momentum, "maxContinuity": continuity,
              "maxAMIMass": ami_mass, "maxGCL": max(gcl, default=math.inf),
              "maxAMIMomentum": ami_momentum,
              "maxAMICoverageAreaError": max(coverage, default=0.0)}
    checks = {
        "reportedConverged": "CONVERGED: all-phase solution written" in log and "FAILED:" not in log,
        "physicalResiduals": all(math.isfinite(x) and x <= 1e-7 for x in (momentum, continuity)),
        "AMIMassConservation": math.isfinite(ami_mass) and ami_mass <= 1e-10,
        "AMIMomentumConservation": math.isfinite(ami_momentum) and ami_momentum <= 1e-10,
        "GCL": len(gcl) == metadata["nPhases"] and all(math.isfinite(x) and x <= 1e-9 for x in gcl),
        "AMICoverage": all(math.isfinite(x) and x <= 1e-6 for x in coverage),
    }
    velocities = [internal_vectors(case / time_name / p["name"] / "U", metadata["nCells"])
                  for p in metadata["phases"]]
    checks["finiteSolution"] = all(math.isfinite(x) for phase in velocities for cell in phase for x in cell)
    pressures = [internal_scalars(case / time_name / p["name"] / "p", metadata["nCells"])
                 for p in metadata["phases"]]
    checks["finitePressure"] = all(math.isfinite(x) for phase in pressures for x in phase)
    if metadata["case"] == "oscillatoryChannel":
        nx, ny = metadata["nx"], metadata["ny"]
        discrete_error = max(abs(v[j*nx+i][0] - phase["discreteUbyY"][j])
                             for v,phase in zip(velocities,metadata["phases"]) for j in range(ny) for i in range(nx))
        continuum_error = max(abs(v[j*nx+i][0] - phase["continuumUbyY"][j])
                              for v,phase in zip(velocities,metadata["phases"]) for j in range(ny) for i in range(nx))
        transverse = max(abs(c[k]) for phase in velocities for c in phase for k in (1,2))
        result.update(discreteMaxAbsError=discrete_error, continuumMaxAbsError=continuum_error,
                      transverseMaxAbsVelocity=transverse)
        checks["independentDiscreteStokesReference"] = discrete_error <= 2e-5 and transverse <= 2e-5
        result["pressureMaxAbs"] = max(abs(x) for phase in pressures for x in phase)
        checks["pressureNullAndCheckerboardSuppression"] = result["pressureMaxAbs"] <= 2e-5
        result["referenceNote"] = "Continuum error includes spatial truncation; pass/fail uses the independent complex-arithmetic finite-volume Stokes reference."
    elif metadata["case"] in ("amiCavity", "amiOscillatoryCavity"):
        phase_difference = max((abs(a-b) for phase in velocities[1:]
                               for c0, c1 in zip(velocities[0],phase) for a,b in zip(c0,c1)), default=0.0)
        result["maxPhaseDifference"] = phase_difference
        if metadata["case"] == "amiCavity":
            checks["steadyPhaseAgreement"] = phase_difference <= 2e-5
        else:
            checks["nonsteadyPhaseResponse"] = phase_difference > 1e-5
        checks["AMIWasPresent"] = len(coverage) == metadata["nPhases"]
        result["referenceNote"] = "Cavity verifies nonlinear convergence, AMI mass conservation and prescribed temporal character; no accuracy benchmark is claimed."
    else:
        raise ValueError("Unknown verification case")
    result["checks"] = checks
    result["passed"] = all(checks.values())
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--time", default="1")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = verify(args.case, args.time)
    text = json.dumps(result, indent=2) + "\n"
    print(text, end="")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text)
    sys.exit(0 if result["passed"] else 1)

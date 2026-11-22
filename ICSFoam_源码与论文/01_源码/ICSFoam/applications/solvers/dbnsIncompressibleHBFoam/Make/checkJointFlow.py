#!/usr/bin/env python3
"""Generate/run a small p/U/q Newton step check; not a propeller validation.
Source env.sh first. --prepare-only generates inputs without invoking OpenFOAM.
The run checks a disturbed uniform-flow state on a 32-cell, three-phase mesh.
"""
import argparse
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def header(name, cls="dictionary"):
    return f"FoamFile {{ version 2.0; format ascii; class {cls}; object {name}; }}\n"


def prepare(case):
    def write(name, content):
        path = case / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    write("system/blockMeshDict", header("blockMeshDict") + """
scale 1;
vertices ((0 0 0)(1 0 0)(1 1 0)(0 1 0)(0 0 1)(1 0 1)(1 1 1)(0 1 1));
blocks (hex (0 1 2 3 4 5 6 7) (4 4 2) simpleGrading (1 1 1));
edges ();
boundary (walls { type wall; faces ((0 4 7 3)(1 2 6 5)(0 1 5 4)(3 7 6 2)(0 3 2 1)(4 5 6 7)); });
mergePatchPairs ();
""")
    write("system/controlDict", header("controlDict") + """
application dbnsIncompressibleHBFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 0.3; deltaT 0.1;
writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 14;
runTimeModifiable false; libs ("libICSFoam.so");
""")
    write("system/fvSchemes", header("fvSchemes") + """
ddtSchemes { default HBM; }
gradSchemes { default Gauss linear; }
divSchemes { default Gauss upwind; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default yes; }
""")
    write("system/fvSolution", header("fvSolution") + """
solvers {}
harmonicBalance { instantsNumber 3; selectedPeriod 1; oversampling false; }
pseudoTime
{
    pseudoTol 1e-10; pseudoTolRel 1e-8; pRefCell 0; pRefValue 0;
    momentumTolerance 1e-10; continuityTolerance 1e-10; fluxConsistencyTolerance 1e-10;
    flowLineSearch true; stageZeroInitialGuess true;
    fullStageCorrectors 1; stageNewtonSteps 3;
    stageKrylovDirections 40; stageKrylovCycles 3; stageKrylovRelTol 0.05;
    stageNewtonRelTol 0.01; stageFiniteDifferenceEpsilon 1e-7;
    stagePreconditioner blockJacobi;
    jointFlowOperator true; jointFlowJacobianAudit true;
    flowFacePseudoTime false; stageAdvectiveMomentum false;
    gaugeReplacesContinuity false; massConsistentMomentum false; pressureFluxCorrectors 0;
    flowVelocityScale 1; flowLengthScale 1; pseudoCoNum 1; pseudoCoNumMax 1;
    flowMaxBacktracks 12;
}
flowSolver { solver GMRES; GMRES { preconditioner LUSGS; nDirections 30; maxIter 1; tolerance 1e-11; relTol 1e-8; } }
""")
    write("constant/HBProperties", header("HBProperties") +
          "wholeMesh { active true; allMesh true; cylCoords false; frequenciesList (6.283185307179586); harmonicsList (1); }\n")
    write("constant/transportProperties", header("transportProperties") +
          "transportModel Newtonian; nu [0 2 -1 0 0 0 0] 0.1;\n")
    write("constant/turbulenceProperties", header("turbulenceProperties") + "simulationType laminar;\n")
    for phase in range(3):
        values = [1e-3 * (1 + phase) * math.sin(2 * math.pi * i / 32) for i in range(32)]
        pressure = (header("p", "volScalarField") +
                    "dimensions [0 2 -2 0 0 0 0];\ninternalField nonuniform List<scalar> 32\n(\n" +
                    "\n".join(map(str, values)) +
                    "\n);\nboundaryField { walls { type zeroGradient; } }\n")
        velocity = (header("U", "volVectorField") +
                    "dimensions [0 1 -1 0 0 0 0]; internalField uniform (1 0.2 0.1);\n" +
                    "boundaryField { walls { type fixedValue; value uniform (1 0.2 0.1); } }\n")
        write(f"0/subTimeLevel{phase}/p", pressure)
        write(f"0/subTimeLevel{phase}/U", velocity)
        if phase == 0:
            write("0/p", pressure)
            write("0/U", velocity)


def check_log(text):
    for required in ("IHBJointNewtonNormCheck", "IHBJointLinearBlocks", "IHBJointTaylor",
                     "IHBJointNewtonAccepted", "IHBFlowAccepted", "IHBFinished converged="):
        if required not in text:
            raise RuntimeError(f"Missing diagnostic: {required}")
    if re.search(r"\b(?:nan|[-+]?inf)\b", text, re.I):
        raise RuntimeError("Non-finite values in solver log")
    if "FOAM FATAL" in text or "IHBStagnated" in text:
        raise RuntimeError("Solver failure or stagnation")
    norm_errors = [float(x) for x in re.findall(r"IHBJointNewtonNormCheck[^\n]*relativeError (\S+)", text)]
    if not norm_errors or max(norm_errors) > 1e-10:
        raise RuntimeError("Packed and stage residual metrics disagree")
    pairs = re.findall(r"IHBFlowAccepted[^\n]* merit (\S+) baselineMerit (\S+)", text)
    if any(float(new) >= float(old) for new, old in pairs):
        raise RuntimeError("Accepted flow step did not reduce the physical merit")
    errors = [float(x) for x in re.findall(r"IHBJointTaylor[^\n]*relativeToLinear (\S+)", text)]
    if len(errors) < 4 or min(errors[:4]) > 1e-2:
        raise RuntimeError("Directional derivative audit failed")
    return {"stepCheckPassed": True, "fullConvergenceRequired": False,
            "acceptedFlowSteps": len(pairs), "maxNormMismatch": max(norm_errors),
            "firstTaylorRelativeErrors": errors[:4]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--binary", type=Path)
    args = parser.parse_args()
    case = Path(tempfile.mkdtemp(prefix="ihb-joint-flow-"))
    prepare(case)
    print(f"Case: {case}", flush=True)
    if args.prepare_only:
        return
    if args.binary is None or not args.binary.is_file():
        parser.error("--binary must name the newly built solver")
    if shutil.which("blockMesh") is None:
        parser.error("Source env.sh: blockMesh is not available")
    for label, command in (("blockMesh", ["blockMesh", "-case", str(case)]),
                           ("solver", [str(args.binary.resolve()), "-case", str(case)])):
        log = case / f"log.{label}"
        with log.open("w") as stream:
            run = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=300)
        allowed = (0, 2) if label == "solver" else (0,)
        if run.returncode not in allowed:
            raise RuntimeError(f"{label} exited {run.returncode}; see {log}")
    result = check_log((case / "log.solver").read_text())
    (case / "check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

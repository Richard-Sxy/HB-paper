#!/usr/bin/env python3
"""Generate fresh, tiny v2112 verification meshes; no external Python packages.

Run blockMesh -case CASE -region phaseN after generation.  HBProperties is
owned by the solver driver and deliberately not generated here.
"""
from pathlib import Path
import argparse
import cmath
import json
import math


def header(kind, name):
    return ("FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            f"    class {kind};\n    object {name};\n}}\n\n")


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def field(name, kind, dimensions, value, patches):
    result = header(kind, name) + f"dimensions {dimensions};\ninternalField uniform {value};\nboundaryField\n{{\n"
    for patch, spec in patches.items():
        result += f"    {patch}\n    {{\n        {spec}\n    }}\n"
    return result + "}\n"


def patch(name, kind, faces, extra=""):
    return (f"    {name}\n    {{\n        type {kind};\n        {extra}\n"
            f"        faces ({' '.join('(' + ' '.join(map(str, f)) + ')' for f in faces)});\n    }}\n")


def mesh_text(vertices, blocks, boundaries):
    vs = "\n".join("    (" + " ".join(map(str, v)) + ")" for v in vertices)
    bs = "\n".join("    hex (" + " ".join(map(str, v)) + ") ("
                   + " ".join(map(str, n)) + ") simpleGrading (1 1 1)" for v, n in blocks)
    return (header("dictionary", "blockMeshDict") + f"scale 1;\nvertices\n(\n{vs}\n);\n"
            f"blocks\n(\n{bs}\n);\nedges ();\nboundary\n(\n{''.join(boundaries)});\nmergePatchPairs ();\n")


def box(x0, x1):
    return [(x0, 0, -0.05), (x1, 0, -0.05), (x1, 1, -0.05), (x0, 1, -0.05),
            (x0, 0, 0.05), (x1, 0, 0.05), (x1, 1, 0.05), (x0, 1, 0.05)]


def channel_mesh(ny):
    return mesh_text(box(0, 1), [(range(8), (2, ny, 1))], [
        patch("left", "cyclic", [(0, 4, 7, 3)], "neighbourPatch right;"),
        patch("right", "cyclic", [(1, 2, 6, 5)], "neighbourPatch left;"),
        patch("walls", "wall", [(0, 1, 5, 4), (3, 7, 6, 2)]),
        patch("frontAndBack", "empty", [(0, 3, 2, 1), (4, 5, 6, 7)]),
    ])


def ami_mesh():
    return mesh_text(box(-1, 0) + box(0, 1), [(range(8), (3, 3, 1)), (range(8, 16), (3, 5, 1))], [
        patch("amiLeft", "cyclicAMI", [(1, 2, 6, 5)],
              "neighbourPatch amiRight; transform noOrdering; AMIMethod faceAreaWeightAMI; lowWeightCorrection -1;"),
        patch("amiRight", "cyclicAMI", [(8, 12, 15, 11)],
              "neighbourPatch amiLeft; transform noOrdering; AMIMethod faceAreaWeightAMI; lowWeightCorrection -1;"),
        patch("lid", "wall", [(3, 7, 6, 2), (11, 15, 14, 10)]),
        patch("walls", "wall", [(0, 4, 7, 3), (0, 1, 5, 4), (8, 9, 13, 12), (9, 10, 14, 13)]),
        patch("frontAndBack", "empty", [(0, 3, 2, 1), (4, 5, 6, 7), (8, 11, 10, 9), (12, 13, 14, 15)]),
    ])


def discrete_channel(ny, viscosity, omega):
    """Solve (i omega I - nu L) uHat = 1 by complex Thomas elimination."""
    a = viscosity * ny * ny
    diagonal = [complex(3 * a if j in (0, ny-1) else 2 * a, omega) for j in range(ny)]
    rhs = [1 + 0j] * ny
    for j in range(1, ny):
        m = -a / diagonal[j-1]
        diagonal[j] -= m * -a
        rhs[j] -= m * rhs[j-1]
    result = [0j] * ny
    result[-1] = rhs[-1] / diagonal[-1]
    for j in range(ny-2, -1, -1):
        result[j] = (rhs[j] + a * result[j+1]) / diagonal[j]
    return result


def generate(destination, nphases=3, ny=8, cases=("oscillatoryChannel", "amiCavity", "amiOscillatoryCavity")):
    viscosity, omega = 0.1, 1.0
    scheme = header("dictionary", "fvSchemes") + """ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; }
"""
    solution = header("dictionary", "fvSolution") + "solvers {}\n"
    control = header("dictionary", "controlDict") + """application fullyImplicitHBFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 1;
deltaT 1;
writeControl timeStep;
writeInterval 1;
writeFormat ascii;
writePrecision 14;
writeCompression off;
timeFormat general;
timePrecision 12;
runTimeModifiable false;
"""
    uhat = discrete_channel(ny, viscosity, omega)
    k = cmath.sqrt(1j * omega / viscosity)
    continuum = [(1 - cmath.cosh(k * ((j+0.5)/ny - 0.5))/cmath.cosh(k/2))/(1j*omega) for j in range(ny)]
    for case_name in cases:
        case = destination / case_name
        is_channel = case_name == "oscillatoryChannel"
        write(case / "system/controlDict", control)
        write(case / "constant/HBProperties", header("dictionary", "HBProperties") + f"""harmonics {(nphases-1)//2};
omega {omega};
nu {viscosity};
velocityScale 1;
lengthScale 1;
rcTime 0.05;
pRefCell 0;
pRefValue 0;
absoluteTolerance 1e-10;
relativeTolerance 1e-9;
physicalTolerance 1e-7;
maxNewton 30;
gmresRestart 300;
gmresMaxIterations 1200;
linearRelativeTolerance 0.0001;
jacobianStep 1e-6;
""")
        write(case / "constant/transportProperties", header("dictionary", "transportProperties") +
              f"transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] {viscosity};\n")
        metadata = {"case": case_name, "nPhases": nphases, "nu": viscosity, "omega": omega,
                    "geometry": "fixed", "nCells": 2*ny if is_channel else 24,
                    "scope": "discrete Stokes HB reference" if is_channel else "nonconformal AMI cavity conservation",
                    "steady": case_name == "amiCavity",
                    "phases": []}
        for phase in range(nphases):
            name = f"phase{phase}"
            theta = 2 * math.pi * phase / nphases
            if is_channel:
                velocity_bcs = {"left": "type cyclic;", "right": "type cyclic;",
                                "walls": "type fixedValue; value uniform (0 0 0);", "frontAndBack": "type empty;"}
                pressure_bcs = {**velocity_bcs, "walls": "type zeroGradient;"}
                force = f"({math.sin(theta):.17g} 0 0)"
                block_mesh = channel_mesh(ny)
                metadata["phases"].append({"name": name, "theta": theta, "bodyForce": math.sin(theta),
                    "discreteUbyY": [(u*cmath.exp(1j*theta)).imag for u in uhat],
                    "continuumUbyY": [(u*cmath.exp(1j*theta)).imag for u in continuum]})
            else:
                lid_speed = 1.0 + (0.2 * math.sin(theta) if case_name == "amiOscillatoryCavity" else 0.0)
                velocity_bcs = {"amiLeft": "type cyclicAMI;", "amiRight": "type cyclicAMI;",
                                "lid": f"type fixedValue; value uniform ({lid_speed:.17g} 0 0);",
                                "walls": "type fixedValue; value uniform (0 0 0);", "frontAndBack": "type empty;"}
                pressure_bcs = {**velocity_bcs, "lid": "type zeroGradient;", "walls": "type zeroGradient;"}
                force = "(0 0 0)"
                block_mesh = ami_mesh()
                metadata["phases"].append({"name": name, "theta": theta, "lidSpeed": lid_speed})
            force_bcs = {p: (spec if "cyclic" in spec or "empty" in spec else "type zeroGradient;")
                         for p, spec in pressure_bcs.items()}
            write(case / f"0/{name}/U", field("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", velocity_bcs))
            write(case / f"0/{name}/p", field("p", "volScalarField", "[0 2 -2 0 0 0 0]", "0", pressure_bcs))
            write(case / f"0/{name}/bodyForce", field("bodyForce", "volVectorField", "[0 1 -2 0 0 0 0]", force, force_bcs))
            write(case / f"system/{name}/blockMeshDict", block_mesh)
            write(case / f"system/{name}/fvSchemes", scheme)
            write(case / f"system/{name}/fvSolution", solution)
        if is_channel:
            metadata.update({"nx": 2, "ny": ny, "height": 1,
                "continuumVsDiscreteMaxError": max(abs(a-b) for a,b in zip(uhat, continuum)),
                "referenceConvention": "U(y,t)=Im[Uhat(y) exp(i omega t)]; forcing=sin(omega t)"})
        write(case / "metadata.json", json.dumps(metadata, indent=2) + "\n")
        allmesh = "#!/bin/sh\nset -eu\ncd \"$(dirname \"$0\")\"\n" + "\n".join(
            f"blockMesh -region phase{i} > log.blockMesh.phase{i} 2>&1" for i in range(nphases)) + "\n"
        write(case / "Allmesh", allmesh)
        (case / "Allmesh").chmod(0o755)
        print(case)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "cases")
    parser.add_argument("--n-phases", type=int, default=3)
    parser.add_argument("--ny", type=int, default=8)
    parser.add_argument("--cases", nargs="+", choices=("oscillatoryChannel", "amiCavity", "amiOscillatoryCavity"),
                        default=["oscillatoryChannel", "amiCavity", "amiOscillatoryCavity"])
    args = parser.parse_args()
    if args.n_phases < 3 or args.n_phases % 2 == 0 or args.ny < 2:
        parser.error("n-phases must be odd and >=3; ny must be >=2")
    generate(args.output.resolve(), args.n_phases, args.ny, args.cases)

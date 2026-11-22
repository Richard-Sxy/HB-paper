#!/usr/bin/env python3
"""Check SST-HB with an independently forced homogeneous periodic solution.

No shear/walls in the analytic cases: F1=0, G=0, and SST reduces to
dk/dt=-betaStar*k*omega+Sk, domega/dt=-beta2*omega**2+Somega.
Forcing uses analytic derivatives, not the implementation's HB matrix.
Wall smoke mode separately exercises native omegaWallFunction/nutLowRe.
All cases are created in a fresh /tmp directory; no real propeller is run.
Source the OpenFOAM v2112 environment and the application's library path first.
"""
import argparse
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def header(name, kind="dictionary"):
    return ("FoamFile {version 2.0; format ascii; class " + kind
            + "; object " + name + ";}\n")


def write_field(path, name, value, vector=False, boundary="zeroGradient", rank=None):
    dim = "[0 1 -1 0 0 0 0]" if vector else {
        "p": "[0 2 -2 0 0 0 0]", "k": "[0 2 -2 0 0 0 0]",
        "omega": "[0 0 -1 0 0 0 0]", "nut": "[0 2 -1 0 0 0 0]"}[name]
    body = f"dimensions {dim};\ninternalField uniform {value};\nboundaryField\n{{\n"
    body += f"walls {{type {boundary};"
    if boundary != "zeroGradient":
        body += f" value uniform {value};"
    body += "}\n"
    if rank is not None:
        body += f"procBoundary{rank}to{1-rank} {{type processor; value uniform {value};}}\n"
    body += "}\n"
    path.mkdir(parents=True, exist_ok=True)
    (path/name).write_text(header(name,"volVectorField" if vector else "volScalarField")+body)


def values(path):
    text=path.read_text()
    match=re.search(r"internalField\s+uniform\s+([^;]+);",text)
    if match:
        return [float(match[1])]
    match=re.search(r"internalField\s+nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)",text,re.S)
    assert match,path
    return list(map(float,match[1].split()))


def setup(root, parallel, perturbed, wall, implicit, rotating):
    template=Path(__file__).resolve().parent/"turbulenceCheck"
    case=root/(("wall" if wall else "perturbed" if perturbed else "exact")
               +("-mpi" if parallel else "-serial")+("-guard" if implicit else "")
               +("-rotating" if rotating else ""))
    case.mkdir()
    (case/"system").mkdir()
    (case/"constant").mkdir()
    common={
        "transportProperties":"transportModel Newtonian; nu [0 2 -1 0 0 0 0] 0.01;",
        "dynamicMeshDict":"dynamicFvMesh staticFvMesh;",
        "turbulenceProperties":("simulationType RAS; RAS {RASModel kOmegaSST; turbulence on; "
                                "printCoeffs on; kOmegaSSTCoeffs {decayControl no;} }"),
        "HBProperties":("wholeMesh {active true; allMesh true; cylCoords false; "
                         "frequenciesList (6.283185307179586); harmonicsList (1);}"),
    }
    if rotating:
        common["dynamicMeshDict"]=(
            "dynamicFvMesh dynamicMotionSolverFvMesh; motionSolverLibs (fvMotionSolvers); "
            "motionSolver solidBody; solidBodyMotionFunction rotatingMotion; "
            "origin (0 0 0); axis (-1 0 0); omega 6.283185307179586;"
        )
    # Native Function1 tables are evaluated under the physical phase clock.
    # The tiny motion probe offset is the existing solver's phase-clock offset.
    probe=1/3/5000
    tables={"k":[],"omega":[]}
    exact=[]
    for j in range(3):
        theta=2*math.pi*j/3
        k=1+0.2*math.sin(theta)
        omega=2+0.3*math.cos(theta)
        exact.append((k,omega))
        tables["k"].append((j/3+probe,0.2*2*math.pi*math.cos(theta)+0.09*k*omega))
        tables["omega"].append((j/3+probe,-0.3*2*math.pi*math.sin(theta)+0.0828*omega**2))
    options="manufactured {type scalarSemiImplicitSource; active true; selectionMode all; volumeMode specific; injectionRateSuSp {"
    for field in ["k","omega"]:
        pairs=[(0,tables[field][0][1])]+tables[field]+[(1+probe,tables[field][0][1])]
        options+=field+" {Su table ("+" ".join(f"({t:.17g} {v:.17g})" for t,v in pairs)+"); Sp constant 0;}"
    common["fvOptions"]=options+"}}"
    if wall:
        # A smoke check only: analytic wall-free forcing is intentionally absent.
        common["fvOptions"]=""
    for name,body in common.items():
        (case/"constant"/name).write_text(header(name)+body+"\n")
    if not parallel:
        shutil.copytree(template/"laminarRegression/constant/polyMesh",case/"constant/polyMesh")
    else:
        for j in range(3):
            for name in common:
                p=case/f"constant/subTimeLevel{j}/{name}"
                p.parent.mkdir(exist_ok=True)
                shutil.copy2(case/"constant"/name,p)
        for rank in range(2):
            shutil.copytree(template/f"perturbed/processor{rank}/constant",case/f"processor{rank}/constant")
    roots=[case] if not parallel else [case/f"processor{rank}" for rank in range(2)]
    for rank,r in enumerate(roots):
        mesh=r/"constant/polyMesh/boundary"
        if not wall:
            mesh.write_text(re.sub(r"\btype\s+wall;","type patch;",mesh.read_text()))
        for j,(k,omega) in enumerate(exact):
            field_dir=r/f"0/subTimeLevel{j}"
            factor=1.02 if perturbed else 1
            n=rank if parallel else None
            write_field(field_dir,"U","(0 0 0)" if wall or rotating else "(1 0 0)",True,"fixedValue",n)
            write_field(field_dir,"p","0",False,"fixedValue",n)
            write_field(field_dir,"k",f"{factor*k:.17g}",False,"fixedValue" if wall else "zeroGradient",n)
            if wall:
                # Match the user's low-Re wall convention, including k=0.
                p=field_dir/"k"
                p.write_text(p.read_text().replace(f"value uniform {factor*k:.17g}","value uniform 0"))
            write_field(field_dir,"omega",f"{factor*omega:.17g}",False,"omegaWallFunction" if wall else "zeroGradient",n)
            write_field(field_dir,"nut",f"{k/omega:.17g}",False,"nutLowReWallFunction" if wall else "calculated",n)
    end=2 if wall or implicit else 2500 if perturbed else 1
    (case/"system/controlDict").write_text(header("controlDict")+f"""
application dbnsIncompressibleHBFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime {end}; deltaT 1;
writeControl timeStep; writeInterval 10000; writeFormat ascii; writePrecision 16;
runTimeModifiable false; libs ("libICSFoam.so");
""")
    method="phaseImplicit" if implicit else "ssprk"
    (case/"system/fvSolution").write_text(header("fvSolution")+f"""
solvers {{ "(k|omega)" {{solver PBiCGStab; preconditioner DILU; tolerance 1e-13; relTol 0;}} }}
harmonicBalance {{instantsNumber 3; selectedPeriod 1; oversampling false;}}
pseudoTime {{
hbTurbulenceMethod {method}; turbulenceHBDamping 1; nHBTurbulenceCorrectors 2;
kTolerance 1e-9; omegaTolerance 1e-9; turbulenceResidualAudit true;
pseudoTol 1e-9; pseudoTolRel 1e-6; pRefCell 0; pRefValue 0;
momentumTolerance 1e-8; continuityTolerance 1e-9; fluxConsistencyTolerance 1e-10;
}}
flowSolver {{solver GMRES; GMRES {{preconditioner LUSGS; nDirections 30; maxIter 1; tolerance 1e-11; relTol 1e-8;}}}}
""")
    (case/"system/fvSchemes").write_text(header("fvSchemes")+"""
ddtSchemes {default HBM;}
gradSchemes {default Gauss linear;}
divSchemes {default Gauss linear; div(phi,U) Gauss upwind; div(phi,k) Gauss upwind; div(phi,omega) Gauss upwind;}
laplacianSchemes {default Gauss linear corrected;}
interpolationSchemes {default linear;}
snGradSchemes {default corrected;}
fluxRequired {default yes;}
wallDist {method meshWave;}
""")
    return case,exact


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--parallel",action="store_true")
    parser.add_argument("--perturbed",action="store_true")
    parser.add_argument("--wall-smoke",action="store_true")
    parser.add_argument("--guard",action="store_true")
    parser.add_argument("--rotating",action="store_true",
                        help="Exercise constant-volume ALE phase geometry; not an AMI test")
    args=parser.parse_args()
    root=Path(tempfile.mkdtemp(prefix="ihb-sst-check-"))
    case,exact=setup(root,args.parallel,args.perturbed,args.wall_smoke,args.guard,args.rotating)
    command=[str(args.binary.resolve()),"-case",str(case)]
    if args.parallel:command=["mpirun","-np","2"]+command+["-parallel"]
    log=case/"log.solver"
    with log.open("w") as stream:
        result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=180)
    text=log.read_text()
    if args.guard:
        assert result.returncode!=0 and "experimental kEpsilon-only path" in text,(result.returncode,log)
        print(json.dumps({"case":str(case),"unsupportedModeRejected":True}))
        return
    if args.wall_smoke:
        assert result.returncode in [0,2] and "IHBSSTResidualAudit" in text and "IHBFinished" in text,(result.returncode,log)
        assert "FOAM FATAL" not in text,log
        print(json.dumps({"case":str(case),"wallSmokeOnly":True,"returncode":result.returncode}))
        return
    assert result.returncode==0 and "IHBFinished converged=1" in text,(result.returncode,log)
    output_roots=[case] if not args.parallel else [case/f"processor{rank}" for rank in range(2)]
    error={"k":0.,"omega":0.}
    for r in output_roots:
        times=[p for p in r.iterdir() if p.is_dir() and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?",p.name)]
        last=max(times,key=lambda p:float(p.name))
        for j,reference in enumerate(exact):
            for field,target in zip(["k","omega"],reference):
                error[field]=max(error[field],max(abs(v-target) for v in values(last/f"subTimeLevel{j}/{field}")))
    assert max(error.values())<2e-7,(error,log)
    report={"case":str(case),"parallel":args.parallel,"perturbed":args.perturbed,
            "rotatingALE":args.rotating,"maxAbsError":error,"physicalConverged":True}
    (case/"check.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()

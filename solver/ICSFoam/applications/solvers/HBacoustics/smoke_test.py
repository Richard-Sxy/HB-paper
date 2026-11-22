"""Run the native incompressible solver on manufactured stationary and rotating panels."""
import math
from pathlib import Path
import shutil
import subprocess
from tempfile import TemporaryDirectory
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[4] / "05_开发"))
from propeller_acoustics.fwh import HarmonicPressure, Medium, Surface, fwh


def foam_file(case, name, body, kind="dictionary"):
    path = case / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"FoamFile {{ version 2.0; format ascii; class {kind}; object {path.name}; }}\n" + body)


def run_case(case, rotating, compressible=False, parallel=False):
    w = 2*math.pi
    speed = w if rotating else 0
    foam_file(case, "system/controlDict", """application dbnsIncompressibleHBFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 0.1;
writeControl timeStep; writeInterval 10; writePrecision 14; libs ("libICSFoam.so");""")
    foam_file(case, "system/fvSolution", """solvers {}
harmonicBalance { instantsNumber 3; selectedPeriod 1; oversampling false; }
pseudoTime { pseudoTol 1e-9; pseudoTolRel 1e-6; pRefCell 0; pRefValue 0;
momentumTolerance 1e-8; continuityTolerance 1e-9; fluxConsistencyTolerance 1e-10; }
flowSolver { solver GMRES; GMRES { preconditioner LUSGS; nDirections 30; maxIter 1;
tolerance 1e-11; relTol 1e-8; } }""")
    foam_file(case, "system/fvSchemes", """ddtSchemes { default HBM; }
gradSchemes { default Gauss linear; } divSchemes { default Gauss linear; div(phi,U) Gauss upwind; }
laplacianSchemes { default Gauss linear corrected; } interpolationSchemes { default linear; }
snGradSchemes { default corrected; } fluxRequired { default yes; }""")
    foam_file(case, "system/blockMeshDict", """scale 1;
vertices ((1 0.2 0.1)(2 0.2 0.1)(2 1.2 0.1)(1 1.2 0.1)
(1 0.2 1.1)(2 0.2 1.1)(2 1.2 1.1)(1 1.2 1.1));
blocks (hex (0 1 2 3 4 5 6 7) (2 1 1) simpleGrading (1 1 1)); edges ();
boundary (emitter { type wall; faces ((1 2 6 5)); }
walls { type wall; faces ((0 4 7 3)(0 1 5 4)(3 7 6 2)(0 3 2 1)(4 5 6 7)); }); mergePatchPairs ();""")
    foam_file(case, "constant/HBProperties", f"wholeMesh {{ active true; allMesh true; "
              f"cylCoords false; frequenciesList ({w}); harmonicsList (1); }}")
    motion = (f"dynamicFvMesh dynamicMotionSolverFvMesh; motionSolverLibs (fvMotionSolvers); "
              f"motionSolver solidBody; solidBodyMotionFunction rotatingMotion; origin (0 0 0); "
              f"axis (0 0 1); omega {w};" if rotating else "dynamicFvMesh staticFvMesh;")
    foam_file(case, "constant/dynamicMeshDict", motion)
    foam_file(case, "constant/transportProperties", "transportModel Newtonian; nu [0 2 -1 0 0 0 0] 1e-6;")
    foam_file(case, "constant/turbulenceProperties", "simulationType laminar;")
    phase_times = np.arange(3)/3
    values = 200+70*np.cos(w*phase_times)+50*np.sin(w*phase_times)
    for i, value in enumerate(values):
        p = value/1000
        foam_file(case, f"0/subTimeLevel{i}/p", f"dimensions [0 2 -2 0 0 0 0]; internalField uniform {p}; "
                  f'boundaryField {{ ".*" {{ type fixedValue; value uniform {p}; }} }}', "volScalarField")
        foam_file(case, f"0/subTimeLevel{i}/U", "dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0); "
                  'boundaryField { ".*" { type fixedValue; value uniform (0 0 0); } }', "volVectorField")
    foam_file(case, "system/hbAcousticsDict", f"""enabled true; wallPatches (emitter);
origin (0 0 0); omega (0 0 {speed}); rhoInf 1000; c0 100; pInf 200;
samplesPerPeriod 360; periods 2; observers {{ receiver {{ position (4 0.7 0.6); }} }}""")
    if rotating:
        with (case / "system/hbAcousticsDict").open("a") as config:
            config.write("\ntones { bladeCount 2; orders 3; referencePressure 1e-6; }\n")
    solver = "dbnsFullyImplicitHBFoam" if compressible else "dbnsIncompressibleHBFoam"
    if compressible:
        # Exercise the termination hook on prescribed wall data. Normalised
        # residuals of the spatially constant state can be unity at roundoff.
        solution = case / "system/fvSolution"
        solution.write_text(solution.read_text().replace("pseudoTol 1e-9;", "pseudoTol 1.01;"))
        tutorial = HERE.parents[2] / "tutorials/shockTube"
        shutil.copyfile(tutorial / "constant/thermophysicalProperties", case / "constant/thermophysicalProperties")
        scheme = case / "system/fvSchemes"
        scheme.write_text(scheme.read_text() + "\nconvectiveFluxScheme { fluxScheme ROE; entropyFixCoeff 0.05; }"
                          "\nviscousFluxScheme { LaxFriedrichJacobian true; }\n")
        for i in range(3):
            foam_file(case, f"0/subTimeLevel{i}/p", "dimensions [1 -1 -2 0 0 0 0]; internalField uniform 101525;"
                      ' boundaryField { ".*" { type fixedValue; value uniform 101525; } }', "volScalarField")
            foam_file(case, f"0/subTimeLevel{i}/T", "dimensions [0 0 0 1 0 0 0]; internalField uniform 300;"
                      ' boundaryField { ".*" { type zeroGradient; } }', "volScalarField")
        acoustic = case / "system/hbAcousticsDict"
        acoustic.write_text(acoustic.read_text().replace("pInf 200;", "pInf 101325;"))
    def run(*args):
        proc = subprocess.run(args, cwd=case, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, timeout=90)
        if proc.returncode or (solver in args and "HB acoustics wrote" not in proc.stdout):
            raise RuntimeError(proc.stdout[-6000:])
        return proc.stdout
    run("blockMesh")
    if parallel:
        foam_file(case, "system/decomposeParDict",
                  "numberOfSubdomains 2; method simple; simpleCoeffs { n (2 1 1); delta 0.001; }")
        for i in range(3):
            for folder in ("constant", "system"):
                shutil.copytree(case / folder, case / folder / f"subTimeLevel{i}",
                                ignore=shutil.ignore_patterns("subTimeLevel*"))
        run("decomposePar")
        for i in range(3):
            run("decomposePar", "-region", f"subTimeLevel{i}")
        command = ["mpirun", "-np", "2", solver, "-parallel"]
    else:
        command = [solver]
    run(*command)
    output, = (case / "acousticData").glob("*-time.dat")
    samples = np.loadtxt(output, skiprows=1)
    assert samples.shape == (720, 2)
    t, measured = samples.T
    np.testing.assert_allclose(t, 1+np.arange(720)/360)
    if rotating:
        surf = Surface([[2, .7, .6]], [[-1, 0, 0]], [0, 0, 0], [0, 0, speed], 0)
        source = HarmonicPressure(phase_times, (values-200)[:, None], 1, 1)
        reference = fwh(surf, source, Medium(1000, 100), [4, .7, .6], t)[:, 2]
        tones, = output.parent.glob("*-bpf.dat")
        numbers = np.loadtxt(tones, skiprows=1, usecols=range(1, 7))
        for order, freq, rms, spl, phase, pref in numbers:
            assert freq == 2*order and pref == 1e-6
            coeff = np.mean(measured*np.exp(-2j*np.pi*freq*(t-t[0])))
            np.testing.assert_allclose(rms, np.sqrt(2)*abs(coeff), rtol=1e-8)
            np.testing.assert_allclose([spl, phase], [20*np.log10(rms/pref), np.angle(coeff, deg=True)], atol=1e-6)
    elif compressible:
        reference = np.full_like(t, -200/(16*np.pi))
    else:
        tau = t-2/100
        pressure = 70*np.cos(w*tau)+50*np.sin(w*tau)
        derivative = w*(-70*np.sin(w*tau)+50*np.cos(w*tau))
        reference = -(derivative/200+pressure/4)/(4*np.pi)
    error = np.linalg.norm(measured-reference)/np.linalg.norm(reference)
    print(f"{solver}, {'rotating' if rotating else 'stationary'} panel, "
          f"{2 if parallel else 1} ranks: {len(t)} samples, relative L2 = {error:.6g}")
    assert error < 1e-3

    if not compressible:
        # Leave time 0 in place to catch accidental reads of an earlier phase field.
        partitions = [case / f"processor{i}" for i in range(2)] if parallel else [case]
        for partition in partitions:
            for k in range(3):
                phase = f"subTimeLevel{k}"
                initial = partition / "0" / phase
                saved = partition / "7" / phase
                saved.mkdir(parents=True)
                shutil.copyfile(initial / "p", saved / "p")
                if rotating:
                    (saved / "polyMesh").mkdir()
                    shutil.copyfile(initial / "polyMesh/points", saved / "polyMesh/points")
        control = case / "system/controlDict"
        control.write_text(control.read_text().replace("startTime 0;", "startTime 7;"))
        settings = case / "system/hbAcousticsDict"
        config = settings.read_text()
        settings.write_text(config + "\noutputName savedPanel;\n")
        log = run(*command, "-replayHBAcoustics")
        assert "IHBAcousticReplay checkpoint 7 phases 3" in log and "Starting pseudotime iteration" not in log
        replay = np.loadtxt(output.parent / "savedPanel-time.dat", skiprows=1)
        np.testing.assert_allclose(replay, samples, atol=1e-10, rtol=1e-10)
        if rotating:
            replay_tones = np.loadtxt(output.parent / "savedPanel-bpf.dat", skiprows=1, usecols=range(1, 7))
            np.testing.assert_allclose(replay_tones, numbers, atol=1e-8, rtol=1e-10)
        def expect_replay_error(message):
            proc = subprocess.run([*command, "-replayHBAcoustics"], cwd=case,
                                  capture_output=True, text=True, timeout=90)
            assert proc.returncode and message in proc.stdout + proc.stderr, proc.stdout + proc.stderr

        if not rotating:
            settings.write_text(config.replace("enabled true;", "enabled false;"))
            expect_replay_error("requires enabled true")
            settings.write_text(config)
        field = "polyMesh/points" if rotating else "p"
        missing = partitions[-1] / "7/subTimeLevel1" / field
        missing.unlink()
        expect_replay_error(str(missing))


if __name__ == "__main__":
    for rotate in (False, True):
        with TemporaryDirectory() as work:
            run_case(Path(work), rotate)
    with TemporaryDirectory() as work:
        run_case(Path(work), False, compressible=True)
    with TemporaryDirectory() as work:
        run_case(Path(work), True, parallel=True)

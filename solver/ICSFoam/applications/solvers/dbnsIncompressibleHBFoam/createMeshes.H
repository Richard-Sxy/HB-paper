// Read the number of instants from fvSolution
fvSolution solutionDict(runTime);

// Read harmonics coeffs dict
dictionary harmonics = solutionDict.subDict("harmonicBalance");

int nO = harmonics.getOrDefault<label>("instantsNumber",3);

PtrList<dynamicFvMesh> subTimeMeshes(nO);

// Resolve filesystem operations against the selected case, including -case runs.
fileName hbCasePath(runTime.path());
hbCasePath.toAbsolute();
const fileName hbConstantPath(hbCasePath/runTime.constant());
const fileName hbSystemPath(hbCasePath/runTime.system());
const fileName hbTempConstantPath(hbCasePath/"tmpCnst");
const fileName hbTempSystemPath(hbCasePath/"tmpSys");

if (!Pstream::parRun())
{
	if(fileHandler().exists(hbTempConstantPath))
	{
		rmDir(hbTempConstantPath);
	}

	if(fileHandler().exists(hbTempSystemPath))
	{
		rmDir(hbTempSystemPath);
	}

	forAll(subTimeMeshes, i)
	{
		fileName subLevConst = hbConstantPath/"subTimeLevel" + Foam::name(i);

		if(fileHandler().exists(subLevConst))
		{
			rmDir(subLevConst);
		}

		fileName subLevSys = hbSystemPath/"subTimeLevel" + Foam::name(i);

		if(fileHandler().exists(subLevSys))
		{
			rmDir(subLevSys);
		}
	}

	if (!cp(hbConstantPath, hbTempConstantPath)
	 || !cp(hbSystemPath, hbTempSystemPath))
	{
		FatalErrorInFunction
			<< "Failed to copy constant/system configuration from case "
			<< hbCasePath << exit(FatalError);
	}

	forAll(subTimeMeshes, i)
	{
		word itrName = Foam::name(i);

		Info<< " \n Create mesh for sub-time level " << itrName << nl << endl;

		if (!cp(hbTempConstantPath, hbConstantPath/"subTimeLevel" + itrName)
		 || !cp(hbTempSystemPath, hbSystemPath/"subTimeLevel" + itrName))
		{
			FatalErrorInFunction
				<< "Failed to copy configuration for sub-time level " << itrName
				<< " in case " << hbCasePath << exit(FatalError);
		}

		subTimeMeshes.set
		(
		  i,
			autoPtr<dynamicFvMesh>
			(
				dynamicFvMesh::New
				(
					args,
					runTime
				)
			)
		);

		subTimeMeshes[i].polyMesh::rename("subTimeLevel" + itrName);
	}
}
else
{
	forAll(subTimeMeshes, i)
	{
		word itrName = Foam::name(i);

		Info<< " \n Create mesh for sub-time level " << itrName << nl << endl;

		subTimeMeshes.set
		(
		  i,
			autoPtr<dynamicFvMesh>
			(
				dynamicFvMesh::New
				(
					args,
					runTime
				)
			)
		);

		subTimeMeshes[i].polyMesh::rename("subTimeLevel" + itrName);
	}
}

IOHBZoneList HB(subTimeMeshes[0]);
HB.setInstants();

scalar Time_temp = runTime.value();

const scalarList& snapshots = HB.selectedSnapshots();

if (nO < 2 || snapshots.size() != nO)
	FatalErrorInFunction << "HB requires instantsNumber >= 2" << exit(FatalError);

const scalar deltaPeriod = snapshots[1] - snapshots[0];
const scalar hbSamplePeriod = nO*deltaPeriod;
// Optional rotor-frame time derivative for blade-passage periods.
const HBRotatingFrame hbFrame(HB, snapshots);
if (hbFrame.active())
{
	// The rotor must advance by exactly 2*pi/Z over the sampled period.
	const scalar turns =
		hbFrame.omega()*hbSamplePeriod/constant::mathematical::twoPi;
	const scalar Z = (mag(turns) > VSMALL) ? 1/mag(turns) : GREAT;
	if (mag(Z - std::round(Z)) > 1e-8*Z || std::round(Z) < 1)
		FatalErrorInFunction << "rotatingFrame: the rotor turns " << turns
			<< " revolutions over the sampled period; it must be 1/Z for an"
			<< " integer Z" << exit(FatalError);
	Info<< "IHBRotatingFrame symmetry order Z " << std::round(Z) << endl;
}

// Guards for the absolute-Cartesian ALE form used by the momentum equation.
const HBZoneList& hbZones = HB;
forAll(hbZones, zonei)
{
	const HBZone& zone = hbZones[zonei];
	const dictionary& zoneDict = HB.subDict(zone.name());
	if (!zoneDict.getOrDefault<bool>("active", true))
		continue;
	// On sliding/moving meshes the stationary cells are unsteady too; a
	// cellZone-only time derivative would silently make them steady.
	if (!zoneDict.getOrDefault<bool>("allMesh", false))
		FatalErrorInFunction << "HB zone " << zone.name()
			<< ": set allMesh true; the incompressible HB time derivative"
			<< " must act on every cell" << exit(FatalError);
	// The momentum source uses Cartesian D*U. ICSFoam cylCoords would
	// differentiate rotating-basis components and drop Omega x U.
	if (zoneDict.get<bool>("cylCoords"))
		FatalErrorInFunction << "HB zone " << zone.name()
			<< ": cylCoords true is not supported by this ALE solver"
			<< exit(FatalError);
	// D is a spectral derivative only if every frequency is periodic over
	// the sampled period; otherwise it is a least-squares operator.
	for (const scalar omega : zone.omegaList())
	{
		const scalar cycles = omega*hbSamplePeriod/constant::mathematical::twoPi;
		if (mag(cycles - std::round(cycles)) > 1e-8*max(scalar(1), mag(cycles)))
			FatalErrorInFunction << "HB zone " << zone.name()
				<< ": frequency " << omega << " rad/s completes " << cycles
				<< " cycles over the sampled period " << hbSamplePeriod
				<< " s; use selectedPeriod = 2*pi*n/omega" << exit(FatalError);
	}
}

forAll(subTimeMeshes, i)
{
	runTime.value()=i*deltaPeriod;
	subTimeMeshes[i].update();
}

// The prescribed motion must repeat after the sampled period, otherwise
// phase N-1 is not followed by phase 0 and D*U mixes inconsistent meshes.
{
	const scalar motionPeriodTolerance = harmonics.getOrDefault<scalar>
		("motionPeriodTolerance", 1e-8);
	fvMesh& mesh0 = subTimeMeshes[0];
	const pointField points0(mesh0.points());
	const label savedIndex = runTime.timeIndex();
	const dimensionedScalar savedDeltaT = runTime.deltaT();
	runTime.setTime(snapshots[0] + hbSamplePeriod, savedIndex + 1);
	runTime.setDeltaT(deltaPeriod, false);
	subTimeMeshes[0].update();
	const boundBox bounds(points0, true);
	const scalar span = max(bounds.mag(), VSMALL);
	// With a rotor frame the rotor returns rotated by 2*pi/Z, the rest unchanged.
	pointField expected(points0);
	if (hbFrame.active())
	{
		const boolList inRotor(hbFrame.pointInRotor(mesh0));
		forAll(expected, pointi)
			expected[pointi] =
				hbFrame.rotate(points0[pointi], hbSamplePeriod, inRotor[pointi]);
	}
	const scalar mismatch =
		points0.size() ? max(mag(mesh0.points() - expected)) : scalar(0);
	const scalar relativeMismatch = returnReduce(mismatch, maxOp<scalar>())/span;
	Info<< "IHBMotionPeriod samplePeriod " << hbSamplePeriod
		<< " relativePointMismatch " << relativeMismatch << endl;
	if (!std::isfinite(relativeMismatch) || relativeMismatch > motionPeriodTolerance)
		FatalErrorInFunction
			<< "Mesh motion is not periodic over the HB sample period "
			<< hbSamplePeriod << " s: relative point mismatch "
			<< relativeMismatch << " > motionPeriodTolerance "
			<< motionPeriodTolerance << ". Match dynamicMeshDict and HB"
			<< " frequencies/selectedPeriod" << exit(FatalError);
	// The phase-motion capture below uses index savedIndex+1 again, so the
	// stored old points of phase 0 remain those at snapshots[0].
	runTime.setTime(Time_temp, savedIndex);
	runTime.setDeltaT(savedDeltaT, false);
}

runTime.value() = Time_temp;

if (!Pstream::parRun())
{
	rmDir(hbTempConstantPath);
	rmDir(hbTempSystemPath);
}


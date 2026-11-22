/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interPhaseChangeFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids with
    phase-change (e.g. cavitation).
    Uses VOF (volume of fluid) phase-fraction based interface capturing.

    The momentum and other fluid properties are of the "mixture" and a
    single momentum equation is solved.

    The set of phase-change models provided are designed to simulate cavitation
    but other mechanisms of phase-change are supported within this solver
    framework.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "phaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"

#include "bubbleCloud.H"
#include "phaseCoupling.H"
#include "ODESolver.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids with"
        " phase-change.\n"
        "Uses VOF (volume of fluid) phase-fraction based interface capturing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"


    phaseCoupling coupling(U,alpha1,p,beta1,bubbles);
    Info<<"creta scaleCouple class successfully"<<endl;
    coupling.initial();

    turbulence->validate();
    volScalarField alpha_1 = alpha1;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    autoPtr<OFstream> RFilePtr_;
    if (Pstream::master())
    {
        RFilePtr_.reset(new OFstream("V.dat"));
    }


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
	
	// Store the particle positions

	// --- Pressure-velocity PIMPLE corrector loop
	while (pimple.loop())
        {
            #include "alphaControls.H"
            volScalarField alpha2_o(alpha2);
            surfaceScalarField rhoPhi
            (
                IOobject
                (
                    "rhoPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimMass/dimTime, Zero)
            );

            mixture->correct();

            #include "alphaEqnSubCycle.H"
            interface.correct();
            
            
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
        bubbles.storeGlobalPositions();


        curlUc = fvc::curl(U);
        
        DUcDt = fvc::ddt(U) + (U & fvc::grad(U));

        GradPc = fvc::grad(p);	
	    // Inject bubbles
        bubbles.move(bubbles,g);

        coupling.update();

        
        alpha1 = min(max(alpha1,scalar(0)),scalar(1));


        scalar volumealpha1 = 0.0;
        forAll(alpha1, celli)
        {
            volumealpha1 += mesh.V()[celli] * (1.0 - alpha1[celli]);
        }
        reduce(volumealpha1, sumOp<scalar>());
        if (Pstream::master())
        {
            RFilePtr_()<<runTime.timeName()<<" "<<volumealpha1<<endl;
        }

	    runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

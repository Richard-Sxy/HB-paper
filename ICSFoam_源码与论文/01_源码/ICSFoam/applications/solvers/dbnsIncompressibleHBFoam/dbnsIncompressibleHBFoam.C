/*---------------------------------------------------------------------------*\

    Copyright (C) 2004-2011 OpenCFD Ltd.
    Copyright (C) 2022 Stefano Oliani

-------------------------------------------------------------------------------
License
    This file is part of ICSFoam.

    ICSFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ICSFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ICSFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "findRefCell.H"
#include <cmath>
#include <cfenv>
#include "sigFpe.H"
#include "pseudotimeControl.H"
#include "dynamicFvMesh.H"
#include "IOMRFCoupledZoneList.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "coupledMatrix.H"
#include "blockFvOperators.H"
#include "LUscalarMatrix.H"

#include "complex.H"
#include "IOHBZoneList.H"
#include "HBRotatingFrame.H"
#include "emptyFvPatch.H"
#include "incompressibleALEFlux.H"
#include "HBKEpsilon.H"
#include "HBKOmegaSST.H"
#include "HBTurbulenceStages.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Fully coupled harmonic-balance solver for incompressible flow.\n"
    );


	#include "addCheckCaseOptions.H"
	#include "setRootCaseLists.H"
	#include "createTime.H"
	#include "createMeshes.H"

	#include "initialise.H"

    bool auditFlowBoundaries=physicalControls.lookupOrDefault<bool>
        ("initialFluxBoundaryAudit",false);
    const label flowRejectedStepLimit=physicalControls.lookupOrDefault<label>
        ("flowRejectedStepLimit",3);
    if (flowRejectedStepLimit<0)
        FatalErrorInFunction << "flowRejectedStepLimit must be nonnegative"
            << exit(FatalError);
    label repeatedMinimumCoRejections=0;


	#include "updateMeshes.H"

    forAll(subTimeMeshes, K)
    {
        turbulence[K].validate();
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting pseudotime iteration loop\n" << endl;
    {
        if (!(flowLineSearch && readInitialAbsoluteFlux))
        {
            #include "refreshPhaseFluxes.H"
        }
        else
        {
            Info<< "IHBRestart preserving saved pressure/velocity/absolute-flux state"
                << endl;
        }
        #include "evaluatePhysicalResiduals.H"
    }

    // Optional diagnostic: evaluate the face-response fixed point without
    // solving or updating cell pressure/velocity. Disabled in production runs.
    const label initialResidualRepeats = physicalControls.lookupOrDefault<label>
        ("initialResidualRepeats", 0);
    if (initialResidualRepeats > 0)
    {
        PtrList<scalarField> savedP(nO), savedNut(nO);
        PtrList<vectorField> savedU(nO);
        forAll(subTimeMeshes,K)
        {
            savedP.set(K,new scalarField(p[K].primitiveField()));
            savedU.set(K,new vectorField(U[K].primitiveField()));
            savedNut.set(K,new scalarField(turbulence[K].nut()().primitiveField()));
        }
        for (label probe=0; probe<initialResidualRepeats; ++probe)
        {
            Info<< "IHBInitialResidualRepeat " << probe+1 << endl;
            #include "evaluatePhysicalResiduals.H"
            scalar dp=0, du=0, dnut=0;
            forAll(subTimeMeshes,K)
            {
                dp=max(dp,gMax(mag(p[K].primitiveField()-savedP[K])));
                du=max(du,gMax(mag(U[K].primitiveField()-savedU[K])));
                dnut=max(dnut,gMax(mag(turbulence[K].nut()().primitiveField()-savedNut[K])));
            }
            Info<< "IHBInitialCellChange p " << dp << " U " << du
                << " nut " << dnut << endl;
        }
    }

    while (runTime.run() && !physicalConverged)
    {
        runTime++;
        flowStepAccepted=true;
        Info<< "Iteration time = " << runTime.timeName() << nl << endl;

        // Retain the original controller bookkeeping and pseudo-Co adaptation.
        // Its linear-residual completion flag does not end the physical solve.
        forAll(subTimeMeshes, K)
        {
            solnControl[K].loop();
        }

        // If only the face-response fixed point is pending, refresh it without
        // asking GMRES to solve an already negligible/zero physical RHS.
        if (!equationsConverged || (flowLineSearch && fluxConsistencyMax>fluxTolerance))
        {
            #include "outerLoop.H"
        }
        if (!flowLineSearch)
        {
            #include "refreshPhaseFluxes.H"
        }
        if (!flowLineSearch || flowStepAccepted)
        {
            #include "evaluatePhysicalResiduals.H"
        }
        else
        {
            // A rejected flow step must not trigger a turbulence/wall update.
            #include "evaluateFlowResiduals.H"
            physicalConverged=physicalConverged && turbulenceConverged;
        }

        if (hbRAS && !freezeTurbulence && !turbulenceConverged && flowStepAccepted)
        {
            for
            (
                label turbulenceCorrector=0;
                turbulenceCorrector<nHBTurbulenceCorrectors && !turbulenceConverged;
                ++turbulenceCorrector
            )
            {
                Info<< "IHBTurbulenceCorrectorStart corrector "
                    << turbulenceCorrector+1 << " total "
                    << nHBTurbulenceCorrectors << endl;
                // Each pass takes a fresh phase snapshot for all three stages.
                #include "correctHBTurbulence.H"
                if (!flowLineSearch)
                {
                    #include "refreshPhaseFluxes.H"
                }
                #include "evaluatePhysicalResiduals.H"
                Info<< "IHBTurbulenceCorrectorEnd corrector "
                    << turbulenceCorrector+1 << " total "
                    << nHBTurbulenceCorrectors << " converged "
                    << turbulenceConverged << endl;
            }
        }

        // U, p and pressure-corrected phi now correspond to the same evaluation.
        runTime.write();
        runTime.printExecutionTime(Info);

        // Repeating the same rejected state at the minimum Co rebuilds the
        // same system. End explicitly instead of spending the iteration budget
        // on identical retries; the final saved state is still un-converged.
        bool allAtMinimumCo=true;
        const scalar minimumCo=physicalControls.lookupOrDefault<scalar>
            ("pseudoCoNumMin",0.1);
        forAll(subTimeMeshes,K)
            allAtMinimumCo=allAtMinimumCo
                && gMax(pseudoCoField[K].primitiveField())<=minimumCo*(1+1e-12);
        repeatedMinimumCoRejections=
            (flowLineSearch && !flowStepAccepted && allAtMinimumCo)
            ? repeatedMinimumCoRejections+1 : 0;
        if (flowRejectedStepLimit>0
         && repeatedMinimumCoRejections>=flowRejectedStepLimit)
        {
            Info<< "IHBStagnated repeated minimum-Co rejected steps "
                << repeatedMinimumCoRejections << "; physical convergence not reached"
                << endl;
            break;
        }

    }

    // Force a final consistent output even if the normal write interval is unmet.
    runTime.writeAndEnd();
    Info<< "IHBFinished converged=" << physicalConverged << endl;
    Info<< "End\n" << endl;

    return physicalConverged ? 0 : 2;
}

// ************************************************************************* //

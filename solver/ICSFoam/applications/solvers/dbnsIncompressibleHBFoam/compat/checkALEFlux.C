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
#include "pseudotimeControl.H"
#include "dynamicFvMesh.H"
#include "IOMRFCoupledZoneList.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "coupledMatrix.H"
#include "blockFvOperators.H"

#include "complex.H"
#include "IOHBZoneList.H"
#include "emptyFvPatch.H"
#include "incompressibleALEFlux.H"
#include "HBKEpsilon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("harmonicCheck", "compare nonzero HB derivative and matrix blocks to analytic sine derivative");
    argList::addBoolOption("wallCheck", "check moving-wall kinematics instead of uniform-flow equilibrium");
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "initialise.H"
    bool auditFlowBoundaries=physicalControls.lookupOrDefault<bool>
        ("initialFluxBoundaryAudit",false);
    #include "updateMeshes.H"
    const bool wallCheck = args.found("wallCheck");
    const bool harmonicCheck=args.found("harmonicCheck");
    const scalar omega=constant::mathematical::twoPi;
    const vector amplitude(0.3,-0.2,0.1);
    if(harmonicCheck)
        forAll(subTimeMeshes,K)
        {
            const vector value=vector(1.2,-0.7,0.4)
                +Foam::sin(omega*(snapshots[K]+hbMotionDeltaT))*amplitude;
            U[K].primitiveFieldRef()=value;
            forAll(U[K].boundaryField(),patchi)
                if(U[K].boundaryField()[patchi].fixesValue())
                    U[K].boundaryFieldRef()[patchi]==value;
        }
    scalar derivativeError=0, blockError=0;

    PtrList<vectorField> initialWall(nO);
    forAll(U,K) initialWall.set(K,new vectorField(U[K].boundaryField()[0]));
    scalar worstPair=0, worstNative=0, worstWall=0, worstMomentum=0, worstMass=0;
    for (label repeat=0;repeat<3;++repeat)
    {
        ++runTime;
        runTime.setDeltaT(0.1*(repeat+1),false);
        // Deliberately trigger v2112's stale native mesh-flux reset.
        forAll(subTimeMeshes,K)
            if(subTimeMeshes[K].moving()) (void)subTimeMeshes[K].phi();
        #include "evaluatePhysicalResiduals.H"
        worstMomentum=max(worstMomentum,momentumMax);
        worstMass=max(worstMass,continuityMax);
        forAll(subTimeMeshes,K)
        {
            const fvMesh& mesh=subTimeMeshes[K];
            if(harmonicCheck)
            {
                const vector derivative=omega*Foam::cos(omega*(snapshots[K]+hbMotionDeltaT))*amplitude;
                vectorField action(mesh.nCells(),Zero);
                forAll(subTimeMeshes,J)
                {
                    blockFvMatrix<vector,tensor> block(mesh);
                    block.diag()=Zero;
                    HB.addBlock(block,K,J);
                    action += (block.diag() & U[J].primitiveField());
                }
                forAll(momSource[K],celli)
                {
                    derivativeError=max(derivativeError,mag(momSource[K][celli]/mesh.V()[celli]+derivative));
                    blockError=max(blockError,mag(action[celli]/mesh.V()[celli]-derivative));
                }
            }
            worstPair=max(worstPair,max(mag(phiAbsolute[K]-phi[K]-meshPhiHB[K])).value());
            if(mesh.moving())
                worstNative=max(worstNative,max(mag(mesh.phi()-meshPhiHB[K])).value());
            if(wallCheck)
            {
                worstWall=max(worstWall,gMax(mag(U[K].boundaryField()[0]-initialWall[K])));
                worstWall=max(worstWall,gMax(mag(phi[K].boundaryField()[0])));
            }
        }
    }
    reduce(derivativeError,maxOp<scalar>()); reduce(blockError,maxOp<scalar>());
    Info<<"ALECheck ranks="<<Pstream::nProcs()<<" wall="<<wallCheck
        <<" pair="<<worstPair<<" native="<<worstNative<<" wallError="<<worstWall
        <<" momentum="<<worstMomentum<<" mass="<<worstMass
        <<" derivativeError="<<derivativeError<<" blockError="<<blockError<<endl;
    bool ok=std::isfinite(worstPair)&&std::isfinite(worstNative)&&std::isfinite(worstWall)
        &&worstPair<1e-12&&worstNative<1e-12&&worstWall<1e-10;
    if(!wallCheck && !harmonicCheck) ok=ok&&std::isfinite(worstMomentum)&&std::isfinite(worstMass)
        &&worstMomentum<1e-8&&worstMass<1e-10;
    ok=ok&&std::isfinite(derivativeError)&&std::isfinite(blockError)
        &&derivativeError<1e-12&&blockError<1e-12;
    (void)physicalConverged;
    return ok?0:2;
}

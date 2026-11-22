#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "coupledMatrix.H"
#include "blockFvOperators.H"
#include "cyclicAMIFvPatch.H"
#include <cmath>
using namespace Foam;
scalar pressureValue(const vector& x) { return Foam::sin(1.3*x.x())+0.3*sqr(x.y())+0.1*x.x()*x.y(); }
vector velocityValue(const vector& x) { return vector(x.x()+0.3*x.y(),0.2*x.x()-x.y(),0.15); }
int main(int argc,char* argv[])
{
    argList::addOption("mode","word","blocking or nonBlocking");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    Pstream::defaultCommsType=args.getOrDefault<word>("mode","nonBlocking")=="blocking"
        ? Pstream::commsTypes::blocking : Pstream::commsTypes::nonBlocking;
    volScalarField p(IOobject("p",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    volVectorField U(IOobject("U",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    mesh.setFluxRequired(p.name());
    scalar worst=0;
    for(label phase=0;phase<3;++phase)
    {
        runTime.setTime(0.037*phase,phase+1);
        runTime.setDeltaT(0.037,false);
        mesh.update();
        forAll(p,i) { p[i]=pressureValue(mesh.C()[i]); U[i]=velocityValue(mesh.C()[i]); }
        forAll(mesh.boundary(),patchi)
        {
            if(p.boundaryField()[patchi].fixesValue())
                forAll(p.boundaryField()[patchi],f)
                    p.boundaryFieldRef()[patchi][f]=pressureValue(mesh.Cf().boundaryField()[patchi][f]);
            if(U.boundaryField()[patchi].fixesValue())
                forAll(U.boundaryField()[patchi],f)
                    U.boundaryFieldRef()[patchi][f]=velocityValue(mesh.Cf().boundaryField()[patchi][f]);
        }
        p.correctBoundaryConditions(); U.correctBoundaryConditions();
        volScalarField rAU(IOobject("rAU",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionedScalar(dimTime,1),zeroGradientFvPatchScalarField::typeName);
        forAll(rAU,i) rAU[i]=1+0.05*mesh.C()[i].x()+0.03*mesh.C()[i].y();
        rAU.correctBoundaryConditions();
        fvScalarMatrix P(-fvm::laplacian(fvc::interpolate(rAU),p));
        const surfaceScalarField pressureFlux(P.flux());
        // Native scalar LDU interface path used by the pressure solver.
        scalarField nativeAction(mesh.nCells(),0);
        P.lduMatrix::Amul(nativeAction,p.primitiveField(),P.boundaryCoeffs(),p.boundaryField().scalarInterfaces(),0);
        const volScalarField nativeDiag(P.A());
        nativeAction += (nativeDiag.primitiveField()*mesh.V()-P.diag())*p.primitiveField();
        fvScalarMatrix insertCopy(P);
        blockFvMatrix<scalar,scalar> S(mesh); S.diag()=Zero; S.source()=Zero;
        S.insertEquation(insertCopy);
        scalarField result(mesh.nCells(),0);
        S.Amul(result,p,mesh);
        scalar ss=gMax(mag(result-nativeAction));
        coupledMatrix joint(mesh,1,1,true);
        joint.dSByS(0,0).diag()=Zero; joint.dSByS(0,0).source()=Zero;
        joint.dSByS(0,0)+=S;
        scalarField assembled(mesh.nCells(),0);
        joint.dSByS(0,0).Amul(assembled,p,mesh);
        scalar assembledError=gMax(mag(assembled-nativeAction));
        const volScalarField fluxDiv(fvc::div(pressureFlux));
        scalar fluxError=gMax(mag(fluxDiv.primitiveField()*mesh.V()-(nativeAction-S.source())));
        scalar coeffError=0, signedFlux=0, absoluteFlux=0;
        label amiFaces=0;
        forAll(mesh.boundary(),patchi)
            if(isA<cyclicAMIFvPatch>(mesh.boundary()[patchi]))
            {
                amiFaces+=mesh.boundary()[patchi].size();
                coeffError=max(coeffError,gMax(mag(S.interfacesUpper()[patchi]+P.boundaryCoeffs()[patchi])));
                signedFlux+=sum(pressureFlux.boundaryField()[patchi]);
                absoluteFlux+=sum(mag(pressureFlux.boundaryField()[patchi]));
            }
        reduce(amiFaces,sumOp<label>()); reduce(coeffError,maxOp<scalar>());
        reduce(signedFlux,sumOp<scalar>()); reduce(absoluteFlux,sumOp<scalar>());
        blockFvMatrix<vector,vector> G(fvm::grad(p)());
        blockFvMatrix<scalar,vector> C(fvm::div(U)());
        vectorField gp(mesh.nCells(),Zero);scalarField cu(mesh.nCells(),0);
        G.Amul(gp,p,mesh);C.Amul(cu,U,mesh);
        const volVectorField gradP(fvc::grad(p));const volScalarField divU(fvc::div(U));
        scalar gradError=gMax(mag(gp-G.source()-mesh.V()*gradP.primitiveField()));
        scalar divError=gMax(mag(cu-C.source()-mesh.V()*divU.primitiveField()));
        scalar scale=max(1.0,gMax(mag(nativeAction)));
        for(const scalar e:{ss,assembledError,fluxError,coeffError,gradError,divError})
        {
            if(!std::isfinite(e)) FatalErrorInFunction<<"Nonfinite audit result"<<exit(FatalError);
            worst=max(worst,e/scale);
        }
        if(!amiFaces) FatalErrorInFunction<<"AMI mesh required"<<exit(FatalError);
        Info<<"AMIBlockAudit phase="<<phase<<" ranks="<<Pstream::nProcs()<<" faces="<<amiFaces
            <<" scalar="<<ss<<" assembled="<<assembledError<<" pressureFlux="<<fluxError
            <<" coefficients="<<coeffError<<" gradient="<<gradError<<" divergence="<<divError
            <<" nativeAMIFluxSum="<<signedFlux<<" nativeAMIFluxL1="<<absoluteFlux<<endl;
    }
    Info<<"AMIBlockAudit worstScaledError="<<worst<<endl;
    return worst<1e-10?0:2;
}

#include "fvCFD.H"
#include "blockFvMatrix.H"
#include "processorFvPatch.H"
#include <cmath>
#include <set>
using namespace Foam;

template<class S, class B>
void allocateInterfaceBlocks(blockFvMatrix<S,B>& A, const fvMesh& mesh)
{
    A.diag() = Zero;
    A.interfacesUpper().resize(mesh.boundary().size());
    A.interfacesLower().resize(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        A.interfacesUpper().set(patchi,new Field<B>(mesh.boundary()[patchi].size(),Zero));
        A.interfacesLower().set(patchi,new Field<B>(mesh.boundary()[patchi].size(),Zero));
    }
}
int main(int argc,char* argv[])
{
    argList::addOption("mode","word","blocking or nonBlocking");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    const word mode=args.getOrDefault<word>("mode","nonBlocking");
    Pstream::defaultCommsType = mode=="blocking"
        ? Pstream::commsTypes::blocking : Pstream::commsTypes::nonBlocking;
    volScalarField q(IOobject("q",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    volVectorField W(IOobject("W",runTime.timeName(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE),mesh);
    blockFvMatrix<scalar,scalar> SS(mesh);
    blockFvMatrix<scalar,vector> SV(mesh);
    blockFvMatrix<vector,vector> VS(mesh);
    blockFvMatrix<vector,tensor> VV(mesh);
    allocateInterfaceBlocks(SS,mesh); allocateInterfaceBlocks(SV,mesh);
    allocateInterfaceBlocks(VS,mesh); allocateInterfaceBlocks(VV,mesh);
    const vector g(1.1,-2.2,0.7);
    const tensor T(1,2,3,4,-1,0.5,0.7,-2,1.2);
    label faces=0, repeatedOwners=0;
    forAll(mesh.boundary(), patchi)
        if(isA<processorFvPatch>(mesh.boundary()[patchi]))
        {
            const labelUList& cells=mesh.boundary()[patchi].faceCells();
            const std::set<label> unique(cells.begin(),cells.end());
            faces+=cells.size(); repeatedOwners+=cells.size()-unique.size();
        }
    reduce(faces,sumOp<label>());
    reduce(repeatedOwners,sumOp<label>());
    if(faces==0) FatalErrorInFunction<<"This check requires processor interfaces"<<exit(FatalError);
    scalar worst=0;
    for(label repeat=0;repeat<3;++repeat)
    {
        forAll(q,celli)
        {
            const vector x=mesh.C()[celli];
            q[celli]=(repeat+1)*(1+2*x.x()+3*x.y()+0.5*x.z());
            W[celli]=vector(x.x()+2*x.y(),3*x.x()-x.y(),2*x.x()+4*x.y()+1)*(repeat+1);
        }
        // Independent oracle: native field evaluation, separate from matrix exchange.
        q.correctBoundaryConditions(); W.correctBoundaryConditions();
        scalarField eSS(mesh.nCells(),0),eSV(mesh.nCells(),0),rSS(mesh.nCells(),0),rSV(mesh.nCells(),0);
        vectorField eVS(mesh.nCells(),Zero),eVV(mesh.nCells(),Zero),rVS(mesh.nCells(),Zero),rVV(mesh.nCells(),Zero);
        forAll(mesh.boundary(), patchi)
        {
            if(!isA<processorFvPatch>(mesh.boundary()[patchi])) continue;
            const scalarField nq(q.boundaryField()[patchi].patchNeighbourField());
            const vectorField nW(W.boundaryField()[patchi].patchNeighbourField());
            const labelUList& cells=mesh.boundary()[patchi].faceCells();
            forAll(cells,facei)
            {
                const vector x=mesh.Cf().boundaryField()[patchi][facei];
                const scalar weight=1+0.13*x.x()+0.21*x.y()+0.07*x.z();
                SS.interfacesUpper()[patchi][facei]=weight;
                SV.interfacesUpper()[patchi][facei]=weight*g;
                VS.interfacesUpper()[patchi][facei]=weight*g;
                VV.interfacesUpper()[patchi][facei]=weight*T;
                eSS[cells[facei]]+=weight*nq[facei];
                eSV[cells[facei]]+=weight*(g & nW[facei]);
                eVS[cells[facei]]+=weight*g*nq[facei];
                eVV[cells[facei]]+=weight*(T & nW[facei]);
            }
        }
        SS.Amul(rSS,q,mesh); SV.Amul(rSV,W,mesh);
        VS.Amul(rVS,q,mesh); VV.Amul(rVV,W,mesh);
        forAll(q,celli)
        {
            for(const scalar error : {mag(rSS[celli]-eSS[celli]),mag(rSV[celli]-eSV[celli]),
                                      mag(rVS[celli]-eVS[celli]),mag(rVV[celli]-eVV[celli])})
            {
                if(!std::isfinite(error)) FatalErrorInFunction<<"Nonfinite result"<<exit(FatalError);
                worst=max(worst,error);
            }
        }
    }
    reduce(worst,maxOp<scalar>());
    Info<<"ProcessorBlockCheck mode="<<mode<<" ranks="<<Pstream::nProcs()
        <<" faces="<<faces<<" repeatedOwners="<<repeatedOwners<<" maxAbsoluteError="<<worst<<endl;
    return worst<=1e-12?0:2;
}

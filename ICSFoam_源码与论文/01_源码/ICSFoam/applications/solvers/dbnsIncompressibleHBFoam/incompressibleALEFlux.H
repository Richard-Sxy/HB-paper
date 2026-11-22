// Shared flux operations for absolute Cartesian velocity on moving HB meshes.
// The phase geometry/mesh flux is fixed during pseudo-time iterations.
#ifndef incompressibleALEFlux_H
#define incompressibleALEFlux_H

namespace Foam
{
// Integrated div(beta*interpolate(deltaU)). Retain the native boundary and
// AMI/processor actions of ICSFoam's velocity divergence, with weighted areas.
inline tmp<blockFvMatrix<scalar,vector>> weightedVelocityDivergence
(
    const volVectorField& U,const surfaceScalarField& beta
)
{
    const fvMesh& mesh=U.mesh();
    const surfaceScalarField& weights=mesh.weights();
    tmp<blockFvMatrix<scalar,vector>> result
        (new blockFvMatrix<scalar,vector>(mesh));
    auto& matrix=result.ref();
    const vectorField area(beta.primitiveField()*mesh.Sf().primitiveField());
    matrix.lower()=-weights.primitiveField()*area;
    matrix.upper()=matrix.lower()+area;
    matrix.negSumDiag();
    matrix.interfacesUpper().resize(mesh.boundary().size());
    matrix.interfacesLower().resize(mesh.boundary().size());
    forAll(U.boundaryField(),patchi)
    {
        const fvPatchVectorField& patch=U.boundaryField()[patchi];
        const auto& w=weights.boundaryField()[patchi];
        const vectorField faceArea
            (beta.boundaryField()[patchi]*mesh.Sf().boundaryField()[patchi]);
        const labelList& cells=patch.patch().faceCells();
        const vectorField lower(-w*faceArea),upper(lower+faceArea);
        matrix.interfacesLower().set(patchi,new vectorField(lower));
        matrix.interfacesUpper().set(patchi,new vectorField(upper));
        const vectorField internal(patch.valueInternalCoeffs(w));
        forAll(cells,i)
            matrix.diag()[cells[i]]+=cmptMultiply(internal[i],faceArea[i]);
        if (!patch.coupled())
        {
            const vectorField boundary(patch.valueBoundaryCoeffs(w));
            forAll(cells,i) matrix.source()[cells[i]]-=boundary[i]&faceArea[i];
        }
    }
    return result;
}

inline void correctAbsoluteFluxBoundary
(
    surfaceScalarField& absolutePhi,
    const volVectorField& U
)
{
    const fvMesh& mesh = U.mesh();
    forAll(mesh.boundary(), patchi)
    {
        // inletOutlet is assignable even though mixed::fixesValue() is true.
        // Preserve its open-boundary flux, following the native HbyA distinction.
        if (!U.boundaryField()[patchi].assignable())
            absolutePhi.boundaryFieldRef()[patchi] =
                U.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
        if (isA<emptyFvPatch>(mesh.boundary()[patchi]))
            absolutePhi.boundaryFieldRef()[patchi] = 0.0;
    }
}

// Read-only diagnostic: pass the absolute flux BEFORE calling the helper.
// All ranks must call this together. It reports physical patches only and
// simulates exactly the current non-assignable/empty override without changing q.
inline void auditAbsoluteFluxBoundary
(
    const surfaceScalarField& absolutePhi,
    const volVectorField& U,
    const volScalarField& p
)
{
    const fvMesh& mesh=U.mesh();
    if (&absolutePhi.mesh()!=&mesh || &p.mesh()!=&mesh)
        FatalErrorInFunction << "Boundary flux audit requires one phase mesh"
            << exit(FatalError);

    forAll(mesh.boundary(),patchi)
    {
        const fvPatch& patch=mesh.boundary()[patchi];
        if (patch.coupled()) continue;
        const fvPatchVectorField& velocityPatch=U.boundaryField()[patchi];
        const fvPatchScalarField& pressurePatch=p.boundaryField()[patchi];
        const scalarField& before=absolutePhi.boundaryField()[patchi];
        scalarField after(before);
        if (!velocityPatch.assignable())
            after=velocityPatch & mesh.Sf().boundaryField()[patchi];
        if (isA<emptyFvPatch>(patch)) after=0.0;

        label faces=before.size();
        scalar beforeSum=0,afterSum=0,sumAbsDifference=0,maxDifferenceByArea=0;
        bool finite=true;
        const scalarField& area=mesh.magSf().boundaryField()[patchi];
        forAll(before,facei)
        {
            const scalar difference=after[facei]-before[facei];
            const scalar speedDifference=mag(difference)/max(area[facei],VSMALL);
            finite=finite && std::isfinite(before[facei])
                && std::isfinite(after[facei]) && std::isfinite(speedDifference);
            beforeSum += before[facei];
            afterSum += after[facei];
            sumAbsDifference += mag(difference);
            maxDifferenceByArea=max(maxDifferenceByArea,speedDifference);
        }
        reduce(faces,sumOp<label>());
        reduce(beforeSum,sumOp<scalar>());
        reduce(afterSum,sumOp<scalar>());
        reduce(sumAbsDifference,sumOp<scalar>());
        reduce(maxDifferenceByArea,maxOp<scalar>());
        reduce(finite,andOp<bool>());
        finite=finite && std::isfinite(beforeSum) && std::isfinite(afterSum)
            && std::isfinite(sumAbsDifference) && std::isfinite(maxDifferenceByArea);
        Info<< "IHBALEBoundaryFluxAudit mesh " << mesh.name()
            << " field " << absolutePhi.name() << " patch " << patch.name()
            << " patchType " << patch.type()
            << " UType " << velocityPatch.type()
            << " pType " << pressurePatch.type()
            << " UFixesValue " << velocityPatch.fixesValue()
            << " UAssignable " << velocityPatch.assignable()
            << " pFixesValue " << pressurePatch.fixesValue()
            << " pAssignable " << pressurePatch.assignable()
            << " faces " << faces << " qBeforeSum " << beforeSum
            << " qAfterSum " << afterSum
            << " overrideSumAbs " << sumAbsDifference
            << " overrideMaxByArea " << maxDifferenceByArea
            << " finite " << finite << endl;
    }
}

// Non-owning list of the HB phase meshes, for coupledMatrix::setEquationMeshes.
// Every coupled solve must apply the AMI of the phase that owns the equation.
template<class MeshType>
inline UPtrList<const fvMesh> hbPhaseMeshes(const PtrList<MeshType>& meshes)
{
    UPtrList<const fvMesh> result(meshes.size());
    forAll(meshes, K) result.set(K, &meshes[K]);
    return result;
}

inline void restoreHBMeshFlux(fvMesh& mesh, const surfaceScalarField& savedFlux)
{
    if (mesh.moving())
    {
        // fvMesh::phi() clears stale flux when the Time index changes.
        // Restore through the non-clearing accessor BEFORE any boundary/model
        // calls (e.g. movingWallVelocity) that may ask for native meshPhi.
        mesh.setPhi() = savedFlux;
        mesh.setPhi().primitiveFieldRef(); // stamp the current pseudo-time index
    }
}
// Keep the phase clock active for the WHOLE assembly/model evaluation:
// patch coefficients may be evaluated lazily inside a discretisation operator.
class HBPhaseTime
{
    Time& time_;
    fvMesh& mesh_;
    const surfaceScalarField& flux_;
    const scalar savedTime_;
    const label savedIndex_;
    const dimensionedScalar savedDeltaT_;
public:
    HBPhaseTime
    (
        Time& runTime, const scalar phaseTime, const scalar motionDeltaT,
        fvMesh& mesh, const surfaceScalarField& savedFlux
    )
    : time_(runTime), mesh_(mesh), flux_(savedFlux),
      savedTime_(runTime.value()), savedIndex_(runTime.timeIndex()),
      savedDeltaT_(runTime.deltaT())
    {
        time_.setTime(phaseTime, savedIndex_);
        time_.setDeltaT(motionDeltaT, false);
        restoreHBMeshFlux(mesh_, flux_);
    }
    ~HBPhaseTime()
    {
        time_.setTime(savedTime_, savedIndex_);
        time_.setDeltaT(savedDeltaT_, false);
        restoreHBMeshFlux(mesh_, flux_);
    }
    HBPhaseTime(const HBPhaseTime&) = delete;
    void operator=(const HBPhaseTime&) = delete;
};
}
#endif

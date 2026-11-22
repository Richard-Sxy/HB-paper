// Optional HB time derivative in the rotor frame for rotor/stator problems that
// repeat after one blade passage T_B=T/Z rather than one revolution T.
//
// U stays the absolute Cartesian velocity. A rotor cell turns with the rotor,
// so its Cartesian U is NOT T_B-periodic, but V=R(t)^T U is (R = rotor rotation).
// With U=R*V:  dU/dt = R dV/dt + Omega x U, and dV/dt|_J = sum_K D_JK V_K, hence
//     dU/dt|_J = sum_K D_JK Q_JK U_K + W U_J,
//     Q_JK = R_J R_K^T = Rot(axis, omega*(t_J-t_K)),  W U = omega*axis x U.
// Stator cells keep the plain D. A rigidly rotating V gives exactly Omega x U
// because the rows of D sum to zero. Scalars (p, k, omega, nut) are frame
// invariant and keep the plain D. Requires a Z-fold symmetric configuration
// (rotor, stator and inflow); this cannot be verified from one mesh alone.
//
// HBProperties, inside the single active HB zone:
//   rotatingFrame { cellZone rotor; origin (0 0 0); axis (0 0 1); omega 125.66; }
#ifndef HBRotatingFrame_H
#define HBRotatingFrame_H

namespace Foam
{

class HBRotatingFrame
{
    bool active_;
    word zoneName_;
    point origin_;
    vector axis_;
    scalar omega_;
    scalarSquareMatrix D_;
    scalarList snapshots_;
    // Per-pair rotor-frame transfer Q_JK and the spin tensor W.
    List<tensor> Q_;
    tensor W_;

    static tensor cross(const vector& a)
    {
        // cross(a) & u == a ^ u
        return tensor(0, -a.z(), a.y(), a.z(), 0, -a.x(), -a.y(), a.x(), 0);
    }

public:

    static tensor rotation(const vector& axis, const scalar angle)
    {
        return Foam::cos(angle)*tensor::I + Foam::sin(angle)*cross(axis)
             + (1 - Foam::cos(angle))*(axis*axis);
    }

    HBRotatingFrame
    (
        const IOHBZoneList& HB, const scalarList& snapshots
    )
    :
        active_(false), origin_(Zero), axis_(Zero), omega_(0),
        snapshots_(snapshots), W_(Zero)
    {
        const HBZoneList& zones = HB;
        label nActive = 0, frameZone = -1;
        forAll(zones, zonei)
        {
            const dictionary& dict = HB.subDict(zones[zonei].name());
            if (!dict.getOrDefault<bool>("active", true)) continue;
            ++nActive;
            if (dict.found("rotatingFrame")) frameZone = zonei;
        }
        if (frameZone < 0) return;
        if (nActive != 1)
            FatalErrorInFunction << "rotatingFrame requires exactly one active"
                << " HB zone" << exit(FatalError);
        const dictionary& dict =
            HB.subDict(zones[frameZone].name()).subDict("rotatingFrame");
        active_ = dict.getOrDefault<bool>("active", true);
        if (!active_) return;
        zoneName_ = dict.get<word>("cellZone");
        origin_ = dict.get<point>("origin");
        axis_ = dict.get<vector>("axis");
        omega_ = dict.get<scalar>("omega");
        if (mag(axis_) <= VSMALL || !std::isfinite(omega_))
            FatalErrorInFunction << "rotatingFrame needs a nonzero axis and"
                << " finite omega" << exit(FatalError);
        axis_ /= mag(axis_);
        D_ = zones[frameZone].D();
        const label nO = snapshots_.size();
        Q_.setSize(nO*nO);
        forAll(snapshots_, J)
            forAll(snapshots_, K)
                Q_[J*nO + K] =
                    rotation(axis_, omega_*(snapshots_[J] - snapshots_[K]));
        W_ = omega_*cross(axis_);
        Info<< "IHBRotatingFrame cellZone " << zoneName_ << " omega " << omega_
            << " axis " << axis_ << " origin " << origin_ << endl;
    }

    bool active() const { return active_; }
    scalar omega() const { return omega_; }

    const tensor& Q(const label J, const label K) const
    {
        return Q_[J*snapshots_.size() + K];
    }

    // Expected position of phase-0 point x after time dt.
    point rotate(const point& x, const scalar dt, const bool inRotor) const
    {
        return inRotor ? origin_ + (rotation(axis_, omega_*dt) & (x - origin_)) : x;
    }

    const labelList& cells(const fvMesh& mesh) const
    {
        const label zoneID = mesh.cellZones().findZoneID(zoneName_);
        if (zoneID < 0)
            FatalErrorInFunction << "rotatingFrame cellZone " << zoneName_
                << " not found" << exit(FatalError);
        return mesh.cellZones()[zoneID];
    }

    boolList pointInRotor(const fvMesh& mesh) const
    {
        boolList inRotor(mesh.nPoints(), false);
        for (const label celli : cells(mesh))
            for (const label pointi : mesh.cellPoints()[celli])
                inRotor[pointi] = true;
        return inRotor;
    }

    // Phase geometry must be the rigid rotation assumed by Q_JK.
    template<class MeshType>
    void checkPhaseGeometry(const PtrList<MeshType>& meshes) const
    {
        if (!active_) return;
        const labelList& rotor = cells(meshes[0]);
        boolList inRotor(meshes[0].nCells(), false);
        for (const label celli : rotor) inRotor[celli] = true;
        scalar mismatch = 0;
        const boundBox bounds(meshes[0].points(), true);
        forAll(meshes, K)
        {
            const vectorField& C0 = meshes[0].C().primitiveField();
            const vectorField& CK = meshes[K].C().primitiveField();
            forAll(CK, celli)
                mismatch = max(mismatch, mag(CK[celli]
                  - rotate(C0[celli], snapshots_[K] - snapshots_[0], inRotor[celli])));
        }
        reduce(mismatch, maxOp<scalar>());
        const scalar relative = mismatch/max(bounds.mag(), VSMALL);
        Info<< "IHBRotatingFrame phaseGeometryMismatch " << relative << endl;
        if (!std::isfinite(relative) || relative > 1e-8)
            FatalErrorInFunction << "Phase meshes are not the rigid rotation"
                << " given by rotatingFrame (relative mismatch " << relative
                << "); match its omega/axis/origin/cellZone to dynamicMeshDict"
                << exit(FatalError);
    }

    // Replace the plain -V*sum_K D_JK U_K of rotor cells by the frame form.
    // Call after HB.addingSource has filled momSource for all cells.
    void correctSource
    (
        PtrList<vectorField>& momSource,
        const UPtrList<volVectorField>& U,
        const fvMesh& mesh
    ) const
    {
        if (!active_) return;
        const scalarField& V = mesh.V();
        const label nO = U.size();
        for (const label celli : cells(mesh))
            for (label J = 0; J < nO; ++J)
            {
                vector rate = W_ & U[J][celli];
                for (label K = 0; K < nO; ++K)
                    rate += D_[J][K]*(Q(J, K) & U[K][celli]);
                momSource[J][celli] = -V[celli]*rate;
            }
    }

    // Rotor-cell HB velocity Jacobian block J,K: V*(D_JK Q_JK + delta_JK W).
    // Call after the plain V*D_JK*I blocks were added for all cells.
    void correctBlock
    (
        tensorField& diag, const label J, const label K, const fvMesh& mesh
    ) const
    {
        if (!active_) return;
        const scalarField& V = mesh.V();
        const tensor frame = D_[J][K]*Q(J, K) + (J == K ? W_ : tensor(Zero));
        const tensor plain = D_[J][K]*tensor::I;
        for (const label celli : cells(mesh))
            diag[celli] += V[celli]*(frame - plain);
    }
};

} // End namespace Foam

#endif

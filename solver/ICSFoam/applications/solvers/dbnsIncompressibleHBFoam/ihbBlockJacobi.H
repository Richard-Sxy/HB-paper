// Application-local block Jacobi for the incompressible all-phase system.
// GPL v3 or later, as ICSFoam.
#ifndef ihbBlockJacobi_H
#define ihbBlockJacobi_H

#include "coupledMatrix.H"
#include "scalarMatrices.H"

namespace Foam
{

class ihbBlockJacobi : public coupledMatrix::preconditioner
{
    const label nScalar_;
    const label nVector_;
    const label nVariables_;
    const scalar pivotTolerance_;
    const bool checkInverse_;
    Field<scalarSquareMatrix> inverseDiagonal_;

    void calculateInverseDiagonal();

public:
    TypeName("ihbBlockJacobi");

    ihbBlockJacobi
    (
        const coupledMatrix::solver& solver,
        const dictionary& controls
    );

    virtual void precondition
    (
        PtrList<scalarField>& scalarResidual,
        PtrList<vectorField>& vectorResidual
    ) const override;
};

}
#endif

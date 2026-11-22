// Application-local block Jacobi for the incompressible all-phase system.
// GPL v3 or later, as ICSFoam.
#include "ihbBlockJacobi.H"
#include <cmath>
#include <limits>
#include "PstreamReduceOps.H"

/**
 * 为多相位压力-速度耦合线性系统实现块 Jacobi 预条件器。
 * 对每个网格单元，收集该单元内各相位压力与速度分量之间的耦合系数，
 * 组装成一个局部稠密矩阵。九相位，每个相位一个压力和三个速度分量时，每个单元的局部矩阵大小为 36*36。
 * 构造阶段对矩阵进行逐行缩放，使用带部分选主元的 LU 分解求出矩阵逆。
 * 应用阶段把局部逆乘到当前单元的残差向量上，再协会各相位的相应分数。
 */

namespace Foam
{

defineTypeNameAndDebug(ihbBlockJacobi, 0);
coupledMatrix::preconditioner::
    adddictionaryConstructorToTable<ihbBlockJacobi>
    addIhbBlockJacobiDictionaryConstructorToTable_;


ihbBlockJacobi::ihbBlockJacobi
(
    const coupledMatrix::solver& solver,
    const dictionary& controls
)
:
    coupledMatrix::preconditioner(solver),
    nScalar_(solver.matrix().nScal()),
    nVector_(solver.matrix().nVect()),
    nVariables_(nScalar_ + 3*nVector_),
    pivotTolerance_
    (
        controls.lookupOrDefault<scalar>("blockJacobiPivotTolerance", 1e-12)
    ),
    checkInverse_(controls.lookupOrDefault<bool>("blockJacobiCheckInverse", false)),
    inverseDiagonal_(solver.matrix().mesh().nCells())
{
    if
    (
        nVariables_ <= 0 || !std::isfinite(pivotTolerance_)
     || pivotTolerance_ <= 0 || pivotTolerance_ >= 1
    )
    {
        FatalErrorInFunction
            << "Require variables and 0 < blockJacobiPivotTolerance < 1"
            << exit(FatalError);
    }
    calculateInverseDiagonal();
}


void ihbBlockJacobi::calculateInverseDiagonal()
{
    const coupledMatrix& matrix = solver_.matrix();
    scalar maximumBackwardError=0;
    scalar maximumScaledIdentityError=0;
    scalar maximumScaledToleranceRatio=0;
    bool inverseCheckFinite=true;
    // Allow roundoff from a small dense solve/product; do not demand exact I.
    const scalar backwardTolerance=
        1000*nVariables_*std::numeric_limits<scalar>::epsilon();

    forAll(inverseDiagonal_, celli)
    {
        scalarSquareMatrix lu(nVariables_, Zero);

        // Flatten [all scalars, xyz(vector 0), xyz(vector 1), ...].
        // Include the diagonal coefficients of ALL variable/phase blocks.
        for (label row=0; row<nScalar_; ++row)
        {
            for (label col=0; col<nScalar_; ++col)
                if (matrix.dSBySExists(row, col))
                    lu[row][col] = matrix.dSByS(row, col).diag()[celli];

            for (label col=0; col<nVector_; ++col)
                if (matrix.dSByVExists(row, col))
                {
                    const vector& block = matrix.dSByV(row, col).diag()[celli];
                    for (direction cmpt=0; cmpt<3; ++cmpt)
                        lu[row][nScalar_+3*col+cmpt] = block.component(cmpt);
                }
        }
        for (label row=0; row<nVector_; ++row)
        {
            for (label col=0; col<nScalar_; ++col)
                if (matrix.dVBySExists(row, col))
                {
                    const vector& block = matrix.dVByS(row, col).diag()[celli];
                    for (direction cmpt=0; cmpt<3; ++cmpt)
                        lu[nScalar_+3*row+cmpt][col] = block.component(cmpt);
                }

            for (label col=0; col<nVector_; ++col)
                if (matrix.dVByVExists(row, col))
                {
                    const tensor& block = matrix.dVByV(row, col).diag()[celli];
                    for (direction i=0; i<3; ++i)
                        for (direction j=0; j<3; ++j)
                            lu[nScalar_+3*row+i][nScalar_+3*col+j] =
                                block.component(3*i+j);
                }
        }

        // Preserve the assembled ORIGINAL coefficients for an independent
        // A*Ainv check. The check does not infer correctness from LU factors.
        autoPtr<scalarSquareMatrix> original;
        if (checkInverse_) original.reset(new scalarSquareMatrix(lu));

        // Scale rows before pivot tests: pressure and velocity rows have
        // different physical units. Invert the original block, not the scaled
        // one: the row scales are also applied to each right-hand side below.
        scalarField rowScale(nVariables_, Zero);
        labelList permutation(nVariables_);
        for (label row=0; row<nVariables_; ++row)
        {
            scalar largest=0;
            permutation[row]=row;
            for (label col=0; col<nVariables_; ++col)
            {
                if (!std::isfinite(lu[row][col]))
                    FatalErrorInFunction
                        << "Nonfinite local block at cell " << celli
                        << ", row " << row << ", column " << col
                        << exit(FatalError);
                largest=max(largest, mag(lu[row][col]));
            }
            if (largest <= VSMALL)
                FatalErrorInFunction
                    << "Zero local block row at cell " << celli
                    << ", row " << row << exit(FatalError);
            rowScale[row]=1.0/largest;
            for (label col=0; col<nVariables_; ++col)
                lu[row][col] *= rowScale[row];
        }

        // Partial-pivot LU of the row-scaled block. Never replace a small
        // pivot by an arbitrary epsilon: reject singular/near-singular blocks.
        for (label col=0; col<nVariables_; ++col)
        {
            label pivot=col;
            scalar largest=mag(lu[col][col]);
            for (label row=col+1; row<nVariables_; ++row)
                if (mag(lu[row][col]) > largest)
                {
                    largest=mag(lu[row][col]);
                    pivot=row;
                }
            if (!std::isfinite(largest) || largest <= pivotTolerance_)
                FatalErrorInFunction
                    << "Singular/near-singular row-scaled block at cell "
                    << celli << ", column " << col << ", pivot " << largest
                    << ", tolerance " << pivotTolerance_
                    << exit(FatalError);
            if (pivot != col)
            {
                for (label j=0; j<nVariables_; ++j)
                {
                    const scalar value=lu[col][j];
                    lu[col][j]=lu[pivot][j];
                    lu[pivot][j]=value;
                }
                const label index=permutation[col];
                permutation[col]=permutation[pivot];
                permutation[pivot]=index;
            }
            for (label row=col+1; row<nVariables_; ++row)
            {
                lu[row][col] /= lu[col][col];
                for (label j=col+1; j<nVariables_; ++j)
                    lu[row][j] -= lu[row][col]*lu[col][j];
            }
        }

        scalarSquareMatrix inverse(nVariables_, Zero);
        scalarField solution(nVariables_, Zero);
        for (label column=0; column<nVariables_; ++column)
        {
            // Solve P*(R*A) x = P*R*e_column through L and U.
            for (label row=0; row<nVariables_; ++row)
            {
                scalar value=permutation[row]==column ? rowScale[column] : 0;
                for (label j=0; j<row; ++j)
                    value -= lu[row][j]*solution[j];
                solution[row]=value;
            }
            for (label row=nVariables_-1; row>=0; --row)
            {
                scalar value=solution[row];
                for (label j=row+1; j<nVariables_; ++j)
                    value -= lu[row][j]*solution[j];
                solution[row]=value/lu[row][row];
                if (!std::isfinite(solution[row]))
                    FatalErrorInFunction
                        << "Nonfinite local inverse at cell " << celli
                        << ", row " << row << ", column " << column
                        << exit(FatalError);
                inverse[row][column]=solution[row];
            }
        }
        inverseDiagonal_[celli]=inverse;

        if (checkInverse_)
        {
            const scalarSquareMatrix& a=original();
            const scalarSquareMatrix& aInverse=inverseDiagonal_[celli];
            scalar aNorm=0, inverseNorm=0, identityError=0;
            scalar scaledANorm=0, scaledInverseNorm=0, scaledIdentityError=0;
            bool productsFinite=true;
            for (label row=0; row<nVariables_; ++row)
            {
                scalar aRow=0, inverseRow=0, errorRow=0;
                scalar scaledARow=0, scaledInverseRow=0, scaledErrorRow=0;
                for (label col=0; col<nVariables_; ++col)
                {
                    aRow += mag(a[row][col]);
                    inverseRow += mag(aInverse[row][col]);
                    scaledARow += mag(rowScale[row]*a[row][col]);
                    // B=R*A, hence B^-1=A^-1*R^-1.
                    scaledInverseRow += mag(aInverse[row][col]/rowScale[col]);
                    scalar product=0, scaledProduct=0;
                    for (label k=0; k<nVariables_; ++k)
                    {
                        product += a[row][k]*aInverse[k][col];
                        scaledProduct += (rowScale[row]*a[row][k])
                            *(aInverse[k][col]/rowScale[col]);
                    }
                    productsFinite=productsFinite && std::isfinite(product)
                        && std::isfinite(scaledProduct);
                    const scalar identity=row==col ? 1 : 0;
                    errorRow += mag(identity-product);
                    scaledErrorRow += mag(identity-scaledProduct);
                }
                aNorm=max(aNorm,aRow);
                inverseNorm=max(inverseNorm,inverseRow);
                identityError=max(identityError,errorRow);
                scaledANorm=max(scaledANorm,scaledARow);
                scaledInverseNorm=max(scaledInverseNorm,scaledInverseRow);
                scaledIdentityError=max(scaledIdentityError,scaledErrorRow);
            }
            const scalar normProduct=aNorm*inverseNorm;
            const scalar scaledNormProduct=scaledANorm*scaledInverseNorm;
            const scalar backwardError=identityError/(normProduct+1);
            // A raw identity-error bound must account for conditioning.
            // Also report the unnormalised scaled error so it remains visible.
            const scalar scaledTolerance=max
            (
                scalar(1e-10), backwardTolerance*(scaledNormProduct+1)
            );
            const scalar scaledRatio=scaledIdentityError/scaledTolerance;
            const bool finite=
                productsFinite && std::isfinite(normProduct) && std::isfinite(scaledNormProduct)
             && std::isfinite(identityError) && std::isfinite(backwardError)
             && std::isfinite(scaledIdentityError) && std::isfinite(scaledRatio);
            inverseCheckFinite=inverseCheckFinite && finite;
            if (finite)
            {
                maximumBackwardError=max(maximumBackwardError,backwardError);
                maximumScaledIdentityError=max
                    (maximumScaledIdentityError,scaledIdentityError);
                maximumScaledToleranceRatio=max
                    (maximumScaledToleranceRatio,scaledRatio);
            }
        }
    }

    if (checkInverse_)
    {
        label checkedCells=inverseDiagonal_.size();
        reduce(checkedCells,sumOp<label>());
        reduce(inverseCheckFinite,andOp<bool>());
        reduce(maximumBackwardError,maxOp<scalar>());
        reduce(maximumScaledIdentityError,maxOp<scalar>());
        reduce(maximumScaledToleranceRatio,maxOp<scalar>());
        // One global summary per constructed preconditioner, never per cell
        // or Krylov application. The option is disabled by default.
        Info<< "IHBBlockJacobiInverseCheck cells " << checkedCells
            << " relativeBackwardErrorMax " << maximumBackwardError
            << " scaledIdentityErrorMax " << maximumScaledIdentityError
            << " scaledToleranceRatioMax " << maximumScaledToleranceRatio
            << " finite " << inverseCheckFinite << endl;
        if
        (
            !inverseCheckFinite || maximumBackwardError>backwardTolerance
         || maximumScaledToleranceRatio>1
        )
            FatalErrorInFunction
                << "Local inverse failed independent A*Ainv verification; "
                << "backward error tolerance " << backwardTolerance
                << ", scaled identity tolerance ratio "
                << maximumScaledToleranceRatio << exit(FatalError);
    }
}


void ihbBlockJacobi::precondition
(
    PtrList<scalarField>& scalarResidual,
    PtrList<vectorField>& vectorResidual
) const
{
    if (scalarResidual.size()!=nScalar_ || vectorResidual.size()!=nVector_)
        FatalErrorInFunction << "Preconditioner variable count mismatch"
            << exit(FatalError);
    forAll(scalarResidual, i)
        if (scalarResidual[i].size()!=inverseDiagonal_.size())
            FatalErrorInFunction << "Scalar preconditioner cell count mismatch"
                << exit(FatalError);
    forAll(vectorResidual, i)
        if (vectorResidual[i].size()!=inverseDiagonal_.size())
            FatalErrorInFunction << "Vector preconditioner cell count mismatch"
                << exit(FatalError);

    scalarField rhs(nVariables_, Zero);
    scalarField result(nVariables_, Zero);
    forAll(inverseDiagonal_, celli)
    {
        for (label i=0; i<nScalar_; ++i) rhs[i]=scalarResidual[i][celli];
        for (label i=0; i<nVector_; ++i)
            for (direction cmpt=0; cmpt<3; ++cmpt)
                rhs[nScalar_+3*i+cmpt]=vectorResidual[i][celli].component(cmpt);

        const scalarSquareMatrix& inverse=inverseDiagonal_[celli];
        for (label row=0; row<nVariables_; ++row)
        {
            scalar value=0;
            for (label col=0; col<nVariables_; ++col)
                value += inverse[row][col]*rhs[col];
            if (!std::isfinite(value))
                FatalErrorInFunction
                    << "Nonfinite preconditioned value at cell " << celli
                    << ", row " << row << exit(FatalError);
            result[row]=value;
        }
        for (label i=0; i<nScalar_; ++i) scalarResidual[i][celli]=result[i];
        for (label i=0; i<nVector_; ++i)
            for (direction cmpt=0; cmpt<3; ++cmpt)
                vectorResidual[i][celli].component(cmpt)=result[nScalar_+3*i+cmpt];
    }
}

}

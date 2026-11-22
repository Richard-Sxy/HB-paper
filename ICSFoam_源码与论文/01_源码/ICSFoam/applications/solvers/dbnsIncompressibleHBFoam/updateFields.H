// Apply the coupled increments directly to the incompressible unknowns.
// p is already kinematic pressure; no density or thermodynamic recovery.
volScalarField& p_ = p[K];
volVectorField& U_ = U[K];

const volScalarField& pIncr = scalarSolIncrement[K];
const volVectorField& UIncr = vectorSolIncrement[K];

forAll(pIncr, celli)
{
    if (!std::isfinite(pIncr[celli]))
        FatalErrorInFunction << "Nonfinite pressure increment" << exit(FatalError);
    for (direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
        if (!std::isfinite(UIncr[celli][cmpt]))
            FatalErrorInFunction << "Nonfinite velocity increment" << exit(FatalError);
}

// Update cell values only. DimensionedField operations retain dimension checks.
// Use a full increment, as in the original solver; pseudo-time controls the step.
p_.ref() += pIncr.internalField();
U_.ref() += UIncr.internalField();

// Boundary values follow their field boundary conditions, not increment values.
p_.correctBoundaryConditions();
U_.correctBoundaryConditions();
MRF[K].correctBoundaryVelocity(U_);

// Do not overwrite the pressure-corrected phi with interpolation of U alone.
// residualsUpdate.H rebuilds phi at the next outer residual evaluation.
// The main loop refreshes phi and physical residuals before stopping/output.

// Shared phase evaluation for assembly and post-update physical checks.
// Spatial momentum block: convection plus the incompressible stress model.
// phi is mesh-relative. Its -meshPhiHB contribution is therefore already
// present in the implicit convection matrix; do not subtract another mesh block.
// Freeze phi and the model coefficients during this linear solve.
// Physical HB and pseudo-time are added separately below.
fvVectorMatrix momentumEqn
(
    fvm::div(phi[K], U_)
  + turbulence[K].divDevReff(U_)
);

// Lagged physical diagonal response for this outer linearisation.
// Exclude pseudo-time so that it does not define the target face flux.
const volScalarField aPhysical(momentumEqn.A());
if (gMin(aPhysical.primitiveField()) <= 0)
{
    FatalErrorInFunction
        << "Positive spatial momentum diagonal required for phase " << K
        << exit(FatalError);
}
volScalarField rAU
(
    IOobject
    (
        "rAU", runTime.timeName(), mesh,
        IOobject::NO_READ, IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rAU", dimTime, 0.0),
    zeroGradientFvPatchScalarField::typeName
);
rAU = 1.0/aPhysical;
rAU.correctBoundaryConditions();
const surfaceScalarField rAUf(fvc::interpolate(rAU));

// Build once: use its flux in the physical residual and its coefficients
// in the compact S block. Pressure-reference treatment is added separately.
fvScalarMatrix pressureBlock(-fvm::laplacian(rAUf, p_));

// Refresh phi and momentumEqn, then evaluate the physical RHS BEFORE
// insertEquation changes matrix sources. No pseudo-time residual is added.
#include "residualsUpdate.H"


// Incompressible physical residuals for phase K (integrated over each cell).
// p is kinematic pressure. No density/energy equation or pseudo-time residual.
// rAU is the spatial response frozen for this outer linearisation.
const volVectorField gradP(fvc::grad(p_));

// Rhie-Chow-type volume flux. pressureBlock is -laplacian(rAUf,p),
// so pressureBlock.flux() supplies the NEGATIVE face pressure-gradient term.
// Using that matrix flux keeps its discretisation consistent with the S block.
if (refreshFaceFlux)
{
    phiAbsolute[K] =
        fvc::flux(U_)
      + fvc::flux(rAU*gradP)
      + pressureBlock.flux();
    correctAbsoluteFluxBoundary(phiAbsolute[K], U_);
    phi[K] = phiAbsolute[K] - meshPhiHB[K];
}

// Evaluate the STORED state. Face-response updates are explicit operations,
// not a side effect of measuring the physical equations.
const volScalarField continuityError(fvc::div(phiAbsolute[K]));

// Rebuild convection from the relative counterpart of that absolute flux.
// The turbulence coefficients and rAU remain frozen in this linearisation.
momentumEqn =
    fvm::div(phi[K], U_)
  + turbulence[K].divDevReff(U_);


// H - A*U includes spatial and boundary sources, divided by cell volume.
// In particular, H() also includes contributions from coupled boundaries.
// Compute this before block insertion mutates the fvMatrix source.
const volVectorField spatialMomentumResidual
(
    momentumEqn.H() - momentumEqn.A()*U_
);

// RHS = minus the physical equation residual.
// momSource[K] already equals -V*sum_j(D_Kj*U_j): do not multiply it by V again.
const scalarField continuityResidual
(
    -mesh.V()*continuityError.primitiveField()
);
const vectorField momentumResidual
(
    mesh.V()
   *(
        spatialMomentumResidual.primitiveField()
      - gradP.primitiveField()
    )
  + momSource[K]
);

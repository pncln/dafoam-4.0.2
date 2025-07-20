// Density equation (continuity equation) for compressible flow

if (printToScreen_)
{
    Info << "Solving rhoEqn" << endl;
}

// Continuity equation for compressible flow
solve
(
    fvm::ddt(rho)
  + fvc::div(phi)
);

// Apply density bounds to prevent negative or unrealistic values
rho = max(rho, dimensionedScalar("rhoMin", rho.dimensions(), 1e-8));
rho = min(rho, dimensionedScalar("rhoMax", rho.dimensions(), 1000.0));

// Update boundary conditions for density
rho.correctBoundaryConditions();

// Store old density for relaxation if needed
volScalarField rhoOld("rhoOld", rho);

// Apply under-relaxation for density in supersonic flows for stability
if (pimple.nCorrPIMPLE() <= 1)
{
    scalar rhoRelaxationFactor = 0.9;
    rho = rhoRelaxationFactor * rho + (1.0 - rhoRelaxationFactor) * rhoOld;
}
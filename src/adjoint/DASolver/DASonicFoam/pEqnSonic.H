volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA("HbyA", U);
HbyA = rAU * UEqn.H();

surfaceScalarField rAUf = fvc::interpolate(rAU);
surfaceScalarField rho_f = fvc::interpolate(rho);
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(rho*HbyA)
  + rho_f * rAUf * fvc::ddtCorr(U, phi)
);

fvScalarMatrix pEqn(
    fvm::ddt(psi, p)
  + fvc::div(phiHbyA)
  - fvm::laplacian(rhorAUf, p)
);

pEqn.solve();

phi = phiHbyA + pEqn.flux();

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Update velocity and kinetic energy
U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

// Update density from equation of state
rho = thermo.rho();

// Bound density if needed
DAUtility::boundVar(allOptions, rho, pimplePrintToScreen);

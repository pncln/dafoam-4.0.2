volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    // Use mass flux for HbyA and proper transient correction to keep dimensions consistent
    fvc::flux(rho*HbyA)
  + ddtPhiCorr(rAU, rho, U, phi)
);

surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi)
    * phiHbyA
);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
        + fvm::div(phid, p)
        - fvm::laplacian(rho*rAU, p)
        ==
        fvOptions(psi, p, rho.name())
    );

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = pEqn.flux();
    }
}

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

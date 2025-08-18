volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(rho*HbyA)
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
        // The momentum equation is formulated in terms of acceleration
        // (pressure gradient divided by density). Consequently the
        // pressure correction uses rAU alone instead of rho*rAU to keep
        // the units consistent
        - fvm::laplacian(rAU, p)
        ==
        fvOptions(psi, p, rho.name())
    );

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        // Combine the predicted flux with the pressure correction flux
        // to maintain a mass flux phi for consistent momentum coupling
        phi = phiHbyA + pEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Update velocity and kinetic energy
// Scale the pressure gradient by rho to obtain acceleration before
// applying the rAU operator
U = HbyA - rAU*(fvc::grad(p)/rho);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

// Update density from equation of state
rho = thermo.rho();

// Bound density if needed
DAUtility::boundVar(allOptions, rho, pimplePrintToScreen);

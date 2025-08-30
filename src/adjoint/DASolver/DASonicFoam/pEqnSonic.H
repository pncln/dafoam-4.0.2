// thermo.rho() may lack mass-density dimensions; reattach them to avoid
// dimension-set inconsistencies in subsequent operations.
rho = thermo.rho()
    * dimensionedScalar(
          "rhoDim",
          dimensionSet(1, -3, 0, 0, 0, 0, 0),
          1.0);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();

surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi)
   * (
        fvc::flux(HbyA)
      + rhorAUf * fvc::ddtCorr(rho, U, phi) / fvc::interpolate(rho)
     )
);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvm::div(phid, p)
      - fvm::laplacian(rhorAUf, p)
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

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

// Keep rho field's physical dimensions consistent after each update
rho = thermo.rho()
    * dimensionedScalar(
          "rhoDim",
          dimensionSet(1, -3, 0, 0, 0, 0, 0),
          1.0);
DAUtility::boundVar(allOptions, rho, pimplePrintToScreen);

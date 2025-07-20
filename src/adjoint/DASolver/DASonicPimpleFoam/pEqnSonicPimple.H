// Pressure equation for compressible flow

if (printToScreen_)
{
    Info << "Solving pEqn" << endl;
}

volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPIMPLE() <= 1)
{
    UEqn.clear();
}

// Calculate face flux for compressible flow
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho)*fvc::flux(HbyA)
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rAU, MRF);

// Cache psi for further use
tmp<volScalarField> tpsi(pThermo.psi());
const volScalarField& psi = tpsi();
tmp<volScalarField> trAU(rAU);

// Compressible pressure equation with acoustic effects
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvm::div(phiHbyA)
      - fvm::laplacian(rho*rAU, p)
     ==
        fvOptions(psi, p, rho.name())
    );

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pEqn.flux();
    }
}

// Explicitly relax pressure for stability in supersonic flow
p.relax();

// Update velocity field
U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// Apply pressure control for closed domains
pressureControl.limit(p);

// Update density based on equation of state
rho = pThermo.rho();

// Apply density bounds
rho = max(rho, dimensionedScalar("rhoMin", rho.dimensions(), 1e-8));

// Calculate and report compressible continuity error
{
    volScalarField contErr = fvc::ddt(rho) + fvc::div(phi);
    scalar sumLocalContErr = contErr.weightedAverage(mesh.V()).value();
    scalar globalContErr = runTime.deltaTValue() * sumLocalContErr;
    cumulativeContErr_ += globalContErr;
    
    if (printToScreen_)
    {
        Info << "time step continuity errors : sum local = " << sumLocalContErr
             << ", global = " << globalContErr
             << ", cumulative = " << cumulativeContErr_
             << endl;
    }
}

// Bound pressure to avoid negative values
p = max(p, pMin);

rho = pThermo.rho();
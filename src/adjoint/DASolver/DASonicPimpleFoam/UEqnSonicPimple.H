// Solve momentum equation for compressible flow with shock capturing

if (printToScreen_)
{
    Info << "Solving UEqn" << endl;
}

// Calculate effective viscosity including artificial viscosity for shock capturing
volScalarField muEff("muEff", turbulence.muEff());

// Add artificial viscosity for shock capturing
if (shockCapturingScheme_ == "artificialViscosity")
{
    volScalarField cellVolume = mesh.V();
    volScalarField cellSize = pow(cellVolume, 1.0/3.0);
    volScalarField artificialViscosity = 
        artificialViscosityCoeff_ * shockSensor * cellSize * rho * mag(U);
    muEff += artificialViscosity;
}

// Momentum equation with compressible effects
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + turbulence.divDevRhoReff(U)
 ==
    fvSource
);

// Add pressure gradient source term for compressible flow
UEqn -= fvc::grad(p);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
              - fvc::snGrad(p)*mesh.magSf()
            )
        )
    );
}

// Store rAU for pressure equation
volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(rho*HbyA));

// Update face fluxes for compressible flow
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Apply under-relaxation for stability in supersonic flow
if (pimple.nCorrPIMPLE() <= 1)
{
    UEqn.relax();
}
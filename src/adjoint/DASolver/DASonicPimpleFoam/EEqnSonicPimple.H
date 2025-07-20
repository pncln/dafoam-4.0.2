// Energy equation for supersonic compressible flow

if (printToScreen_)
{
    Info << "Solving EEqn" << endl;
}

// Get references to energy field
volScalarField& he = pThermo.he();

// Calculate effective thermal diffusivity including shock capturing
volScalarField alphaEff("alphaEff", turbulence.alphaEff());

// Add artificial thermal diffusivity for shock capturing
if (shockCapturingScheme_ == "artificialViscosity")
{
    volScalarField cellVolume = mesh.V();
    volScalarField cellSize = pow(cellVolume, 1.0/3.0);
    volScalarField artificialAlpha = 
        artificialViscosityCoeff_ * shockSensor * cellSize * rho * mag(U) * pThermo.Cp();
    alphaEff += artificialAlpha;
}

// Energy equation including kinetic energy effects for supersonic flow
fvScalarMatrix EEqn
(
    fvm::ddt(rho, he) + fvm::div(phi, he)
  + fvc::ddt(rho, 0.5*magSqr(U)) + fvc::div(phi, 0.5*magSqr(U))
  + (
        he.name() == "e"
      ? fvc::div
        (
            fvc::absolute(phi/fvc::interpolate(rho), U),
            p,
            "div(phiv,p)"
        )
      : -dpdt
    )
  - fvm::laplacian(alphaEff, he)
 ==
    rho*(U&g)
  + rad.Sh(pThermo, he)
  + fvSourceEnergy
  + fvOptions(rho, he)
);

EEqn.relax();

fvOptions.constrain(EEqn);

EEqn.solve();

fvOptions.correct(he);

// Update thermodynamic properties
pThermo.correct();

// Apply temperature bounds for stability in supersonic flow
const volScalarField& T = pThermo.T();
volScalarField& TRef = const_cast<volScalarField&>(T);
TRef = max(TRef, dimensionedScalar("TMin", T.dimensions(), 50.0));
TRef = min(TRef, dimensionedScalar("TMax", T.dimensions(), 5000.0));
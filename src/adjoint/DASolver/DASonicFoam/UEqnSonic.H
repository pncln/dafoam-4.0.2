// Solve the Momentum equation

// MRF.correctBoundaryVelocity(U);   // MRF is not used
if (hasFvSource_)
{
    daFvSourcePtr_->calcFvSource(fvSource);
}

tmp<fvVectorMatrix> tUEqn(
    // Build the conventional momentum equation. The body-force fvSource is
    // defined per unit volume and will be added explicitly below.
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + daTurbulenceModelPtr_->divDevRhoReff(U));
fvVectorMatrix& UEqn = tUEqn.ref();

// Constrain by fvOptions if available
fvOptions.constrain(UEqn);

// Add external body-force source term
if (hasFvSource_)
{
    UEqn -= fvSource;
}

UEqn.relax();

// Solve the momentum equation and retrieve solver performance
SolverPerformance<vector> solverU = solve(UEqn == -fvc::grad(p));

// Apply any fvOptions corrections to the solved field
fvOptions.correct(U);

DAUtility::primalResidualControl(solverU, pimplePrintToScreen, "U", daGlobalVarPtr_->primalMaxRes);

// bound U
DAUtility::boundVar(allOptions, U, pimplePrintToScreen);

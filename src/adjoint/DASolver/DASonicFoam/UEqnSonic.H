// Solve the Momentum equation

// MRF.correctBoundaryVelocity(U);   // MRF is not used
if (hasFvSource_)
{
    daFvSourcePtr_->calcFvSource(fvSource);
}

// Build the conventional momentum equation. The body-force fvSource will be
// added explicitly via the source term below.
tmp<fvVectorMatrix> tUEqn(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + daTurbulenceModelPtr_->divDevRhoReff(U));
fvVectorMatrix& UEqn = tUEqn.ref();

// Add external body-force source term explicitly to the equation source
if (hasFvSource_)
{
    UEqn.source() -= fvSource;
}

UEqn.relax();

// Apply any fvOptions constraints after relaxation
fvOptions.constrain(UEqn);

// Solve the momentum equation and retrieve solver performance
SolverPerformance<vector> solverU = solve(UEqn == -fvc::grad(p));

// Apply any fvOptions corrections to the solved field
fvOptions.correct(U);

DAUtility::primalResidualControl(solverU, pimplePrintToScreen, "U", daGlobalVarPtr_->primalMaxRes);

// bound U
DAUtility::boundVar(allOptions, U, pimplePrintToScreen);

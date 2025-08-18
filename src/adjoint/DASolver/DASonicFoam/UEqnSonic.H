// Solve the Momentum equation

// MRF.correctBoundaryVelocity(U);   // MRF is not used 
if (hasFvSource_)
{
    daFvSourcePtr_->calcFvSource(fvSource);
}


tmp<fvVectorMatrix> tUEqn(
    // Build the conventional momentum equation. The body-force fvSource is
    // supplied in acceleration units (force per unit mass), so convert it to
    // a force per unit volume by multiplying with rho before adding it to
    // the equation.
    fvm::ddt(rho, U)
    + fvm::div(phi, U)
    + daTurbulenceModelPtr_->divDevRhoReff(U)
    - rho*fvSource);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

// get the solver performance info such as initial
// and final residuals
// Solve the momentum equation
SolverPerformance<vector> solverU = solve(UEqn == -fvc::grad(p));

DAUtility::primalResidualControl(solverU, pimplePrintToScreen, "U", daGlobalVarPtr_->primalMaxRes);
 
								 

// bound U
DAUtility::boundVar(allOptions, U, pimplePrintToScreen);

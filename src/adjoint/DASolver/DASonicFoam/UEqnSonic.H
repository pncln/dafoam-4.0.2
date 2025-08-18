// Solve the Momentum equation

// MRF.correctBoundaryVelocity(U);   // MRF is not used 
if (hasFvSource_)
{
    daFvSourcePtr_->calcFvSource(fvSource);
}


tmp<fvVectorMatrix> tUEqn(
    // Divide the conventional momentum equation by rho so that all terms
    // have the dimensions of acceleration. This allows the body-force
    // fvSource, which is supplied in acceleration units, to be added
    // directly without any additional scaling.
    (fvm::ddt(rho, U)
    + fvm::div(phi, U)
    + daTurbulenceModelPtr_->divDevRhoReff(U))/rho
    // fvSource already represents an acceleration (force per unit mass)
    // and thus can be combined with the above terms directly
    - fvSource);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

// get the solver performance info such as initial
// and final residuals
// Divide the pressure gradient by rho so that the RHS has
// the dimensions of acceleration, matching the LHS
// of the momentum equation in DASonicFoam
SolverPerformance<vector> solverU = solve(UEqn == -fvc::grad(p)/rho);

DAUtility::primalResidualControl(solverU, pimplePrintToScreen, "U", daGlobalVarPtr_->primalMaxRes);
 
								 

// bound U
DAUtility::boundVar(allOptions, U, pimplePrintToScreen);

// Solve the Momentum equation

// MRF.correctBoundaryVelocity(U);   // MRF is not used 
if (hasFvSource_)
{
    daFvSourcePtr_->calcFvSource(fvSource);
}
fvVectorMatrix UEqn(
    fvm::ddt(rho, U)
    + fvm::div(phi, U)
    + turbulencePtr_->divDevRhoReff(U)
 ==
    fvSource);
				
UEqn.relax();

// get the solver performance info such as initial
// and final residuals
SolverPerformance<vector> solverU = solve(UEqn == -fvc::grad(p));

DAUtility::primalResidualControl(solverU, pimplePrintToScreen, "U", daGlobalVarPtr_->primalMaxRes);
 
								 

// bound U
DAUtility::boundVar(allOptions, U, pimplePrintToScreen);

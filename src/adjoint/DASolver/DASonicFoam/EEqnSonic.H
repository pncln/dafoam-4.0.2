{
    // Calculate energy source term if actuator disk is present
    if (hasFvSource_)
    {
        fvSourceEnergy = fvSource & U;
    }

    // Sync thermo's T with our T field before energy equation
    if (&(thermo.T()) != &T) { thermo.T() = T; }

    // Access thermo.he() directly
    volScalarField& e = thermo.he();

    fvScalarMatrix EEqn(
        fvm::ddt(rho, e)
      + fvm::div(phi, e)
      + fvc::ddt(rho, K)
      + fvc::div(phi, K)
      + (e.name() == "e"
          ? fvc::div(fvc::absolute(phi / fvc::interpolate(rho), U), p, "div(phiv,p)")
          : -dpdt)
      - fvm::laplacian(turbulencePtr_->alphaEff(), e)
      - fvSourceEnergy);
EEqn.relax();

    SolverPerformance<scalar> solverE = EEqn.solve();

    DAUtility::primalResidualControl(solverE, pimplePrintToScreen, "e", daGlobalVarPtr_->primalMaxRes);

    // Bound internal energy
    DAUtility::boundVar(allOptions, e, pimplePrintToScreen);

    // Update thermodynamic properties
    thermo.correct();
    
    // After thermo correction, update our T field to match thermo's T
    if (&T != &(thermo.T())) { T = thermo.T(); }
}

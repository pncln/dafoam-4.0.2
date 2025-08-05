{
    // Calculate energy source term if actuator disk is present
    if (hasFvSource_)
    {
        fvSourceEnergy = fvSource & U;
    }

    // Sync thermo fields with our state variables before energy equation
    thermo.p() = p;
    thermo.T() = T;

    // Access thermo.he() directly
    volScalarField& e = thermo.he();

    fvScalarMatrix EEqn(
        fvm::ddt(rho, e) + fvm::div(phi, e)
        + fvc::ddt(rho, K) + fvc::div(phi, K)
        + fvc::div(fvc::absolute(phi/fvc::interpolate(rho), U), p, "div(phiv,p)")
        - fvm::laplacian(turbulence.alphaEff(), e)
        ==
        fvOptions(rho, e)
        + fvSourceEnergy
    );

    EEqn.relax();

    SolverPerformance<scalar> solverE = EEqn.solve();

    DAUtility::primalResidualControl(solverE, pimplePrintToScreen, "e", daGlobalVarPtr_->primalMaxRes);

    // Bound internal energy
    DAUtility::boundVar(allOptions, e, pimplePrintToScreen);

    // Update thermodynamic properties
    thermo.correct();
}

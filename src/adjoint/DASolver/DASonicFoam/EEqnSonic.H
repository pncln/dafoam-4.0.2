{


    fvScalarMatrix EEqn(
       
        fvm::ddt(rho, e) + fvm::div(phi, e)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + fvc::div(fvc::absolute(phi/fvc::interpolate(rho), U), p, "div(phiv,p)")
	  - fvm::laplacian(turbulence.alphaEff(), e)

     ==
        fvOptions(rho, e)
    );

    EEqn.relax();

    // get the solver performance info such as initial
    // and final residuals
    SolverPerformance<scalar> solverE = EEqn.solve();

    DAUtility::primalResidualControl(solverE, pimplePrintToScreen, "e", daGlobalVarPtr_->primalMaxRes);

    // bound e
    DAUtility::boundVar(allOptions, e, pimplePrintToScreen);

    thermo.correct();
}

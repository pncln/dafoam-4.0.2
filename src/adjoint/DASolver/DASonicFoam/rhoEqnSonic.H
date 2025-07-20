/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    This file is modified from OpenFOAM's source code
    src/finiteVolume/cfdTools/compressible/rhoEqn.H

    OpenFOAM: The Open Source CFD Toolbox

    Copyright (C): 2011-2016 OpenFOAM Foundation

    OpenFOAM License:

        OpenFOAM is free software: you can redistribute it and/or modify it
        under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.
    
        OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
        ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
        FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
        for more details.
    
        You should have received a copy of the GNU General Public License
        along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

{
    // Density continuity equation - exactly like original sonicFoam
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
      - fvc::div(fvSource)
    );

    // get the solver performance info such as initial
    // and final residuals
    SolverPerformance<scalar> solverRho = rhoEqn.solve();
    
    DAUtility::primalResidualControl(solverRho, pimplePrintToScreen, "rho", daGlobalVarPtr_->primalMaxRes);
    
    // bound rho
    DAUtility::boundVar(allOptions, rho, pimplePrintToScreen);
}

// ************************************************************************* //
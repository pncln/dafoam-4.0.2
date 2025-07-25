/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for DASonicFoam
    
    This class is modified from OpenFOAM's source code
    applications/solvers/compressible/sonicFoam

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

#ifndef DASonicFoam_H
#define DASonicFoam_H

#include "DASolver.H"
#include "addToRunTimeSelectionTable.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControlDF.H"
#include "DAFvSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DASonicFoam Declaration
\*---------------------------------------------------------------------------*/

class DASonicFoam
    : public DASolver
{

protected:
    /// pimple pointer
    autoPtr<pimpleControlDF> pimplePtr_;

    /// thermal model
    autoPtr<psiThermo> pThermoPtr_;

    /// pressure field pointer
    autoPtr<volScalarField> pPtr_;

    /// density field pointer
    autoPtr<volScalarField> rhoPtr_;

    /// velocity field pointer
    autoPtr<volVectorField> UPtr_;

    /// surface flux field pointer
    autoPtr<surfaceScalarField> phiPtr_;

    /// pressure time deriv pointer
    autoPtr<volScalarField> dpdtPtr_;

    /// kinetic energy pointer
    autoPtr<volScalarField> KPtr_;

    /// turbulence pointer
    autoPtr<compressible::turbulenceModel> turbulencePtr_;

    /// DATurbulenceModel pointer
    autoPtr<DATurbulenceModel> daTurbulenceModelPtr_;

    /// DASource pointer
    autoPtr<DAFvSource> daFvSourcePtr_;

    /// fvSource term
    autoPtr<volVectorField> fvSourcePtr_;

    /// fvSource term for the energy equation
    autoPtr<volScalarField> fvSourceEnergyPtr_;

    /// whether to have fvSource term
    label hasFvSource_ = 0;

    /// continuity error
    scalar cumulativeContErr_ = 0.0;

    /// whether to write mesh for the reduceIO
    label reduceIOWriteMesh_ = 0;

public:
    TypeName("DASonicFoam");
    // Constructors

    //- Construct from components
    DASonicFoam(
        char* argsAll,
        PyObject* pyOptions);

    //- Destructor
    virtual ~DASonicFoam()
    {
    }

    /// initialize fields and variables
    virtual void initSolver();

    /// solve the primal equations
    virtual label solvePrimal();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
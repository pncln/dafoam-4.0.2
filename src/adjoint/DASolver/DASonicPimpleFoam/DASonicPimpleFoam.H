/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for DASonicPimpleFoam - supersonic compressible flow solver
        with shock capturing capabilities
    
    This class is modified from DAPimpleFoam to handle compressible supersonic flows
    Based on OpenFOAM's sonicFoam and rhoPimpleFoam solvers

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

#ifndef DASonicPimpleFoam_H
#define DASonicPimpleFoam_H

#include "DASolver.H"
#include "addToRunTimeSelectionTable.H"
#include "fluidThermo.H"
#include "compressibleTurbulenceModel.H"
#include "pimpleControlDF.H"
#include "DAFvSource.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DASonicPimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DASonicPimpleFoam
    : public DASolver
{

protected:
    /// pimple pointer
    autoPtr<pimpleControlDF> pimplePtr_;

    /// pressure field pointer
    autoPtr<volScalarField> pPtr_;

    /// velocity field pointer
    autoPtr<volVectorField> UPtr_;

    /// density field pointer
    autoPtr<volScalarField> rhoPtr_;

    /// surface flux field pointer
    autoPtr<surfaceScalarField> phiPtr_;

    /// thermo pointer
    autoPtr<fluidThermo> pThermoPtr_;

    /// turbulence pointer
    autoPtr<compressible::turbulenceModel> turbulencePtr_;

    /// DATurbulenceModel pointer
    autoPtr<DATurbulenceModel> daTurbulenceModelPtr_;

    /// DASource pointer
    autoPtr<DAFvSource> daFvSourcePtr_;

    /// fvSource term
    autoPtr<volVectorField> fvSourcePtr_;

    /// whether to have fvSource term
    label hasFvSource_ = 0;

    /// continuity error
    scalar cumulativeContErr_ = 0.0;

    /// pressure control for closed domains
    autoPtr<pressureControl> pressureControlPtr_;

    /// the primal IO interval
    label IOInterval_ = 1;

    /// whether to write mesh for the reduceIO
    label reduceIOWriteMesh_ = 0;

    /// energy source term pointer
    autoPtr<volScalarField> fvSourceEnergyPtr_;

    /// whether to have energy source term
    label hasEnergySource_ = 0;

    /// shock capturing scheme name
    word shockCapturingScheme_;

    /// shock sensor field
    autoPtr<volScalarField> shockSensorPtr_;

    /// artificial viscosity coefficient
    scalar artificialViscosityCoeff_;

    /// whether to use local time stepping
    label useLocalTimeStep_;

public:
    TypeName("DASonicPimpleFoam");
    // Constructors

    //- Construct from components
    DASonicPimpleFoam(
        char* argsAll,
        PyObject* pyOptions);

    //- Destructor
    virtual ~DASonicPimpleFoam()
    {
    }

    /// initialize fields and variables
    virtual void initSolver();

    /// solve the primal equations
    virtual label solvePrimal();

    /// solve the adjoint equation using the fixed-point iteration method
    virtual label solveAdjointFP(
        Vec dFdW,
        Vec psi);

    /// calculate adjoint residuals
    scalar calcAdjointResiduals(
        const double* psi,
        const double* dFdW,
        double* adjRes);

protected:
    /// calculate shock sensor for artificial viscosity
    void calcShockSensor();

    /// apply artificial viscosity for shock capturing
    void applyArtificialViscosity();

    /// bound density and pressure to physical values
    void boundThermodynamicProperties();

    /// calculate time step for local time stepping
    void calcLocalTimeStep();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
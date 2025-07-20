/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for DASonicPimpleFoam - compressible supersonic flow solver

\*---------------------------------------------------------------------------*/

#ifndef DAResidualSonicPimpleFoam_H
#define DAResidualSonicPimpleFoam_H

#include "DAResidual.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"
#include "adjustPhi.H"
#include "fluidThermo.H"
#include "compressibleTurbulenceModel.H"
#include "pressureControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAResidualSonicPimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DAResidualSonicPimpleFoam
    : public DAResidual
{
protected:
    /// \name These are state variables, state residuals, and partial derivatives
    //@{
    volVectorField& U_;
    volVectorField URes_;

    volScalarField& p_;
    volScalarField pRes_;

    volScalarField& rho_;
    volScalarField rhoRes_;

    volScalarField& T_;
    volScalarField TRes_;

    surfaceScalarField& phi_;
    surfaceScalarField phiRes_;
    //@}

    /// fvSource term for momentum
    volVectorField& fvSource_;

    /// fvSource term for energy
    volScalarField& fvSourceEnergy_;

    /// DATurbulenceModel object
    DATurbulenceModel& daTurb_;

    /// thermodynamics object
    fluidThermo& thermo_;

    /// pimpleControl object which will be initialized in this class
    pimpleControl pimple_;

    /// pressure control for closed domains
    pressureControl& pressureControl_;

    /// whether to has fvSource term
    label hasFvSource_ = 0;

    /// whether to has energy source term
    label hasEnergySource_ = 0;

    /// shock sensor field
    volScalarField& shockSensor_;

    /// artificial viscosity coefficient
    scalar artificialViscosityCoeff_;

    /// shock capturing scheme
    word shockCapturingScheme_;

public:
    TypeName("DASonicPimpleFoam");
    // Constructors

    //- Construct from components
    DAResidualSonicPimpleFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    //- Destructor
    virtual ~DAResidualSonicPimpleFoam()
    {
    }

    // Members

    /// clear the members
    virtual void clear();

    /// compute residual
    virtual void calcResiduals(const dictionary& options);

    /// update any intermediate variables that are dependent on state variables and are used in calcResiduals
    virtual void updateIntermediateVariables();

    /// update the boundary condition for all the states in the selected solver
    virtual void correctBoundaryConditions();

    /// calculate preconditioner matrix for PETSc
    virtual void calcPCMatWithFvMatrix(Mat PCMat);

protected:
    /// calculate shock sensor for artificial viscosity
    void calcShockSensor();

    /// apply artificial viscosity for shock capturing
    void applyArtificialViscosity();

    /// bound thermodynamic properties to physical values
    void boundThermodynamicProperties();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
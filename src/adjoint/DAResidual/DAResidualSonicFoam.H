/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for DASonicFoam

\*---------------------------------------------------------------------------*/

#ifndef DAResidualSonicFoam_H
#define DAResidualSonicFoam_H

#include "DAResidual.H"
#include "addToRunTimeSelectionTable.H"
#include "pimpleControl.H"
#include "adjustPhi.H"
#include "thermodynamicConstants.H"
#include "pressureControl.H"
#include "constrainPressure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAResidualSonicFoam Declaration
\*---------------------------------------------------------------------------*/

class DAResidualSonicFoam
    : public DAResidual
{

protected:
    /// \name These are state variables, state residuals, and partial derivatives
    //@{
    volVectorField& U_;
    volVectorField URes_;

    volScalarField& p_;
    volScalarField pRes_;

    volScalarField& e_;
    volScalarField ERes_;

    surfaceScalarField& phi_;
    surfaceScalarField phiRes_;

    volScalarField rhoRes_;
    //@}

    /// fvSource term
    volVectorField& fvSource_;

    /// fvSource term for the energy equation
    volScalarField& fvSourceEnergy_;

    /// thermophysical property
    psiThermo& thermo_;

    /// \name These are intermediate variables
    //@{
    volScalarField& rho_;
    volScalarField& alphat_;
    volScalarField& psi_;
    volScalarField& dpdt_;
    volScalarField& K_;
    //@}

    /// DATurbulenceModel object
    DATurbulenceModel& daTurb_;

    /// pimpleControl object which will be initialized in this class
    pimpleControl pimple_;

    /// \name These are constants to update the intermediate variables
    //@{
    scalar molWeight_;
    scalar Cp_;
    //@}

    /// whether to have fvSource term
    label hasFvSource_ = 0;

public:
    TypeName("DASonicFoam");
    // Constructors

    //- Construct from components
    DAResidualSonicFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    //- Destructor
    virtual ~DAResidualSonicFoam()
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

    virtual void calcPCMatWithFvMatrix(Mat PCMat);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
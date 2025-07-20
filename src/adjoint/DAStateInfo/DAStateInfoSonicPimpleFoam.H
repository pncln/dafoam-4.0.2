/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for DASonicPimpleFoam - compressible supersonic flow solver

\*---------------------------------------------------------------------------*/

#ifndef DAStateInfoSonicPimpleFoam_H
#define DAStateInfoSonicPimpleFoam_H

#include "DAStateInfo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAStateInfoSonicPimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class DAStateInfoSonicPimpleFoam
    : public DAStateInfo
{

public:
    TypeName("DASonicPimpleFoam");
    // Constructors

    //- Construct from components
    DAStateInfoSonicPimpleFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAStateInfoSonicPimpleFoam()
    {
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
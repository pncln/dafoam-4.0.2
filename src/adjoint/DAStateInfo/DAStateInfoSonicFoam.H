/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for DASonicFoam

\*---------------------------------------------------------------------------*/

#ifndef DAStateInfoSonicFoam_H
#define DAStateInfoSonicFoam_H

#include "DAStateInfo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAStateInfoSonicFoam Declaration
\*---------------------------------------------------------------------------*/

class DAStateInfoSonicFoam
    : public DAStateInfo
{

public:
    TypeName("DASonicFoam");
    // Constructors

    //- Construct from components
    DAStateInfoSonicFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAStateInfoSonicFoam()
    {
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
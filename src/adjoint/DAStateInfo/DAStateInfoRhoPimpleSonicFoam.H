/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for DARhoPimpleSonicFoam

\*---------------------------------------------------------------------------*/

#ifndef DAStateInfoRhoPimpleSonicFoam_H
#define DAStateInfoRhoPimpleSonicFoam_H

#include "DAStateInfo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAStateInfoRhoPimpleSonicFoam Declaration
\*---------------------------------------------------------------------------*/

class DAStateInfoRhoPimpleSonicFoam
    : public DAStateInfo
{

public:
    TypeName("DARhoPimpleSonicFoam");
    // Constructors

    //- Construct from components
    DAStateInfoRhoPimpleSonicFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAStateInfoRhoPimpleSonicFoam()
    {
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
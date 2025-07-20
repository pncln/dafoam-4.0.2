/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAStateInfoRhoPimpleSonicFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAStateInfoRhoPimpleSonicFoam, 0);
addToRunTimeSelectionTable(DAStateInfo, DAStateInfoRhoPimpleSonicFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAStateInfoRhoPimpleSonicFoam::DAStateInfoRhoPimpleSonicFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel)
    : DAStateInfo(modelType, mesh, daOption, daModel)
{
    stateInfo_["volScalarStates"].append("p");
    stateInfo_["volScalarStates"].append("he");  
    stateInfo_["modelStates"].append("nut");
    stateInfo_["volVectorStates"].append("U");
    stateInfo_["surfaceScalarStates"].append("phi");

    // correct the names for model states based on the selected physical model at runtime
    daModel.correctModelStates(stateInfo_["modelStates"]);

    /* 
    Description:
        Adjoint state connectivity info for supersonic flow
    
                     U      he     p     nut    phi
        URes         2      2      1      1      0
        ERes         2      2      2      1      0
        pRes         3      2      2      2      1
        phiRes       2      2      1      1      0
    */

    stateResConInfo_.set(
        "URes",
        {
            {"U", "p", "he", "nut", "phi"}, // lv0
            {"U", "p", "he", "nut"}, // lv1
            {"U", "he"} // lv2
        });

    stateResConInfo_.set(
        "ERes",
        {
            {"U", "p", "he", "nut", "phi"}, // lv0
            {"U", "p", "he", "nut"}, // lv1
            {"U", "p", "he"} // lv2
        });

    stateResConInfo_.set(
        "pRes",
        {
            {"U", "p", "he", "nut", "phi"}, // lv0
            {"U", "p", "he", "nut", "phi"}, // lv1
            {"U", "p", "he", "nut"}, // lv2
            {"U"} // lv3
        });

    stateResConInfo_.set(
        "phiRes",
        {
            {"U", "p", "he", "nut", "phi"}, // lv0
            {"U", "p", "he", "nut"}, // lv1
            {"U", "he"}, // lv2
        });

    // need to correct connectivity for physical models for each residual
    daModel.correctStateResidualModelCon(stateResConInfo_["URes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["ERes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["pRes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["phiRes"]);

    // add physical model residual connectivity
    daModel.addModelResidualCon(stateResConInfo_);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
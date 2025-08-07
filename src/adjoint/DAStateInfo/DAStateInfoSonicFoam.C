/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAStateInfoSonicFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAStateInfoSonicFoam, 0);
addToRunTimeSelectionTable(DAStateInfo, DAStateInfoSonicFoam, dictionary);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

DAStateInfoSonicFoam::DAStateInfoSonicFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel)
    : DAStateInfo(modelType, mesh, daOption, daModel)
{
    /*
    Description:
        Register the names of state variables
        NOTE:
        For model variables, such as turbulence model, register specific names
        For example, register "nut" to modelStates for RANS turbulence models,
        Then, we will call correctModelStates(stateInfo_["modelStates"]) to modify
        "nut" based on the selected turbulence model. For example, for SA model,
        correctModelStates will just replace "nut" with "nuTilda", for SST model,
        it will replace "nut" with "k" and append "omega" to modelStates.
        In other words, the model variables will be modified based on the selected
        models at runtime.
    */

    // States (harmonized with DARhoPimpleFoam)
    stateInfo_["volScalarStates"].append("p");
    stateInfo_["volScalarStates"].append("T");
    stateInfo_["modelStates"].append("nut");
    stateInfo_["volVectorStates"].append("U");
    stateInfo_["surfaceScalarStates"].append("phi");

    // correct model-state names (e.g., "nut" -> "nuTilda" for SA)
    daModel.correctModelStates(stateInfo_["modelStates"]);

    /* 
    Description:
        Adjoint state connectivity info, numbers denote the level of connectivity
        N/A means this state does not connect to the corrsponding residual 
    
                     U      T      p     nut    phi
        URes         2      2      1      1      0
        TRes         2      2      2      1      0
        pRes         3      2      2      2      1
        phiRes       2      2      1      1      0
    */

    stateResConInfo_.set(
        "URes",
        {
            {"U", "p", "T", "nut", "phi"}, // lv0
            {"U", "p", "T", "nut"}, // lv1
            {"U", "T"} // lv2
        });

    stateResConInfo_.set(
        "TRes",
        {
            {"U", "p", "T", "nut", "phi"}, // lv0
            {"U", "p", "T", "nut"}, // lv1
            {"U", "p", "T"} // lv2
        });

    stateResConInfo_.set(
        "pRes",
        {
            {"U", "p", "T", "nut", "phi"}, // lv0
            {"U", "p", "T", "nut", "phi"}, // lv1
            {"U", "p", "T", "nut"}, // lv2
            {"U"} // lv3
        });

    stateResConInfo_.set(
        "phiRes",
        {
            {"U", "p", "T", "nut", "phi"}, // lv0
            {"U", "p", "T", "nut"}, // lv1
            {"U", "T"}, // lv2
        });

    // correct connectivity for physical models for each residual
    daModel.correctStateResidualModelCon(stateResConInfo_["URes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["TRes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["pRes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["phiRes"]);

    // add physical model residual connectivity
    daModel.addModelResidualCon(stateResConInfo_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

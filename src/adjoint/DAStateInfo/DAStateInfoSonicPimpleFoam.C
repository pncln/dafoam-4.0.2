/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAStateInfoSonicPimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAStateInfoSonicPimpleFoam, 0);
addToRunTimeSelectionTable(DAStateInfo, DAStateInfoSonicPimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAStateInfoSonicPimpleFoam::DAStateInfoSonicPimpleFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel)
    : DAStateInfo(modelType, mesh, daOption, daModel)
{
    /*
    Description:
        Register the names of state variables for compressible supersonic flow
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

    // Compressible flow state variables
    stateInfo_["volScalarStates"].append("p");      // pressure
    stateInfo_["volScalarStates"].append("rho");    // density
    stateInfo_["volScalarStates"].append("T");      // temperature
    stateInfo_["modelStates"].append("nut");        // turbulent viscosity
    stateInfo_["volVectorStates"].append("U");      // velocity
    stateInfo_["surfaceScalarStates"].append("phi"); // face flux

    // correct the names for model states based on the selected physical model at runtime
    daModel.correctModelStates(stateInfo_["modelStates"]);

    /* 
    Description:
        Adjoint state connectivity info for compressible supersonic flow
        Numbers denote the level of connectivity
        N/A means this state does not connect to the corresponding residual 
    
                     U      T      p     rho    nut    phi
        URes         2      2      1      1      1      0
        TRes         2      2      2      2      1      0
        pRes         3      2      2      2      2      1
        rhoRes       2      2      1      2      1      1
        phiRes       2      2      1      1      1      0
    
        ******************************** NOTE 1 **********************************
        One does not need to specify connectivity for each physical model, set the 
        connectivity for original variables instead. For example, for turbulence models,
        set nut. Then, how is nut connected to the other turbulence states will be 
        set in the DAModel class. This is done by calling correctStateResidualModelCon. 
        For example, for SA model we just replace nut with nuTilda, for SST model, we need 
        to add extract connectivity since nut depends on grad(U), k, and omega. We need
        to do this for other physical models such as radiation models.
        **************************************************************************
    
        ******************************** NOTE 2 **********************************
        Do not specify physical model connectivity here, because they will be added
        by calling addModelResidualCon. For example, for the SA turbulence
        model, it will add the nuTildaRes to stateResConInfo_ and setup
        its connectivity automatically.
        **************************************************************************

    */

    // Momentum equation connectivity
    stateResConInfo_.set(
        "URes",
        {
            {"U", "p", "rho", "T", "nut", "phi"}, // lv0
            {"U", "p", "rho", "T", "nut"}, // lv1
            {"U", "rho", "T"} // lv2
        });

    // Pressure equation connectivity
    stateResConInfo_.set(
        "pRes",
        {
            {"U", "p", "rho", "T", "nut", "phi"}, // lv0
            {"U", "p", "rho", "T", "nut", "phi"}, // lv1
            {"U", "p", "rho", "T", "nut"}, // lv2
            {"U", "rho", "T"} // lv3
        });

    // Density equation connectivity
    stateResConInfo_.set(
        "rhoRes",
        {
            {"U", "p", "rho", "T", "nut", "phi"}, // lv0
            {"U", "p", "rho", "T", "nut"}, // lv1
            {"rho", "T"}, // lv2
        });

    // Face flux connectivity
    stateResConInfo_.set(
        "phiRes",
        {
            {"U", "p", "rho", "T", "nut", "phi"}, // lv0
            {"U", "p", "rho", "T", "nut"}, // lv1
            {"U", "rho", "T"}, // lv2
        });

    // Temperature/Energy equation connectivity
    stateResConInfo_.set(
        "TRes",
        {
            {"U", "T", "p", "rho", "nut", "phi"}, // lv0
            {"U", "T", "p", "rho", "nut"}, // lv1
            {"T", "rho"} // lv2
        });

    // need to correct connectivity for physical models for each residual
    daModel.correctStateResidualModelCon(stateResConInfo_["URes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["pRes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["rhoRes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["phiRes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["TRes"]);

    // add physical model residual connectivity
    daModel.addModelResidualCon(stateResConInfo_);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
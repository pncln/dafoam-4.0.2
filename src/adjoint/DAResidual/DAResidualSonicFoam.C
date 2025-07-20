/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAResidualSonicFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAResidualSonicFoam, 0);
addToRunTimeSelectionTable(DAResidual, DAResidualSonicFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAResidualSonicFoam::DAResidualSonicFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAResidual(modelType, mesh, daOption, daModel, daIndex),
      // initialize and register state variables and their residuals, we use macros defined in macroFunctions.H
      setResidualClassMemberVector(U, dimensionSet(1, -2, -2, 0, 0, 0, 0)),
      setResidualClassMemberScalar(p, dimensionSet(1, -3, -1, 0, 0, 0, 0)),
      setResidualClassMemberScalar(e, dimensionSet(0, 2, -2, 0, 0, 0, 0)),
      setResidualClassMemberPhi(phi),
      rhoRes_(
          IOobject(
              "rhoRes",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          mesh,
          dimensionedScalar("rhoRes", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0),
          "zeroGradient"),
      fvSource_(const_cast<volVectorField&>(
          mesh_.thisDb().lookupObject<volVectorField>("fvSource"))),
      fvSourceEnergy_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("fvSourceEnergy"))),
      thermo_(const_cast<psiThermo&>(
          mesh_.thisDb().lookupObject<psiThermo>("thermophysicalProperties"))),
      rho_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("rho"))),
      alphat_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("alphat"))),
      psi_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("thermo:psi"))),
      dpdt_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("dpdt"))),
      K_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("K"))),
      daTurb_(const_cast<DATurbulenceModel&>(daModel.getDATurbulenceModel())),
      // create pimpleControl
      pimple_(const_cast<fvMesh&>(mesh), "PIMPLE")
{

    // initialize fvSource
    const dictionary& allOptions = daOption.getAllOptions();
    if (allOptions.subDict("fvSource").toc().size() != 0)
    {
        hasFvSource_ = 1;
    }

    // get molWeight and Cp from thermophysicalProperties
    const IOdictionary& thermoDict = mesh.thisDb().lookupObject<IOdictionary>("thermophysicalProperties");
    dictionary mixSubDict = thermoDict.subDict("mixture");
    dictionary specieSubDict = mixSubDict.subDict("specie");
    molWeight_ = specieSubDict.getScalar("molWeight");
    dictionary thermodynamicsSubDict = mixSubDict.subDict("thermodynamics");
    Cp_ = thermodynamicsSubDict.getScalar("Cp");

    if (daOption_.getOption<label>("debug"))
    {
        Info << "molWeight " << molWeight_ << endl;
        Info << "Cp " << Cp_ << endl;
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DAResidualSonicFoam::clear()
{
    /*
    Description:
        Clear all members to avoid memory leak because we will initalize 
        multiple objects of DAResidual. Here we need to delete all members
        in the parent and child classes
    */
    URes_.clear();
    pRes_.clear();
    eRes_.clear();
    phiRes_.clear();
}

void DAResidualSonicFoam::calcResiduals(const dictionary& options)
{
    /*
    Description:
        This is the function to compute residuals.
    
    Input:
        options.isPC: 1 means computing residuals for preconditioner matrix.
        This essentially use the first order scheme for div(phi,U), div(phi,e)

        p_, he_, U_, phi_, etc: State variables in OpenFOAM
    
    Output:
        URes_, pRes_, eRes_, phiRes_: residuals
    */

    // We dont support transonic = false
    label transonic = 0;

    // ********* U Residuals **********
    // copied and modified from UEqn.H
    fvVectorMatrix UEqn(
        fvm::ddt(rho_, U_)
        + fvm::div(phi_, U_)
        + daTurb_.divDevRhoReff(U_)
        - fvSource_);

    URes_ = UEqn & U_;

    // ********* E Residuals **********
    // Energy equation exactly like original sonicFoam
    fvSourceEnergy_ = fvSource_ & U_;

    fvScalarMatrix EEqn(
        fvm::ddt(rho_, e_) + fvm::div(phi_, e_)
        + fvc::ddt(rho_, K_) + fvc::div(phi_, K_)
        + fvc::div(fvc::absolute(phi_/fvc::interpolate(rho_), U_), p_, "div(phiv,p)")
        - fvm::laplacian(daTurb_.alphaEff(), e_)
        - fvSourceEnergy_);

    eRes_ = EEqn & e_;

    // ********* rho Residuals **********
    // Density continuity equation exactly like original sonicFoam
    fvScalarMatrix rhoEqn(
        fvm::ddt(rho_)
        + fvc::div(phi_)
        - fvSource_);

    rhoRes_ = rhoEqn & rho_;

    // ********* p Residuals **********
    // copied and modified from pEqn.H
    rho_ = thermo_.rho();
    volScalarField rAU(1.0 / UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_ * rAU));
    
    volVectorField HbyA("HbyA", U_);
    HbyA = rAU * UEqn.H();

    surfaceScalarField phid(
        "phid",
        fvc::interpolate(psi_)
       *(
            fvc::flux(HbyA)
          + rhorAUf*fvc::ddtCorr(rho_, U_, phi_)/fvc::interpolate(rho_)
        )
    );

    // we do not support transonic
    if (transonic)
    {
        FatalErrorIn("") << "We do not support transonic = true" << abort(FatalError);
    }

    fvScalarMatrix pEqn(
        fvm::ddt(psi_, p_)
        + fvm::div(phid, p_)
        - fvm::laplacian(rhorAUf, p_));

    pRes_ = pEqn & p_;

    // ********* phi Residuals **********
    phiRes_ = pEqn.flux() - phi_;

    // need to normalize residuals
    normalizeResiduals(URes_);
    normalizeResiduals(pRes_);
    normalizeResiduals(eRes_);
    normalizeResiduals(rhoRes_);
    normalizeResiduals(phiRes_);
}

void DAResidualSonicFoam::updateIntermediateVariables()
{
    /*
    Description:
        Update the intermediate variables that depend on the state variables.
        This function will be called when the state variables are updated.
        So we need to update intermediate variables such as rho, psi, etc.
    */

    // update the density based on the equation of state
    rho_ = thermo_.rho();

    // update psi
    psi_ = thermo_.psi();

    // update alphat
    daTurb_.correctAlphat();

    // update thermodynamic properties
    thermo_.correct();

    // update kinetic energy
    K_ = 0.5 * magSqr(U_);

    // update pressure time derivative
    if (thermo_.dpdt())
    {
        dpdt_ = fvc::ddt(p_);
    }
}

void DAResidualSonicFoam::correctBoundaryConditions()
{
    /*
    Description:
        Update the boundary condition for all the states in the selected solver
    */

    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    he_.correctBoundaryConditions();
    // phi is surface field, no need to correct boundary conditions
}

void DAResidualSonicFoam::calcPCMatWithFvMatrix(Mat PCMat)
{
    /*
    Description:
        Calculate the preconditioner matrix using fvMatrix
    */

    // create pimpleControl
    pimpleControl pimple(const_cast<fvMesh&>(mesh_));

    dictionary pcOptions;
    pcOptions.set("isPC", 1);
    this->calcResiduals(pcOptions);

    // URes
    this->setResidualToPCMat(URes_, "U", PCMat);

    // pRes
    this->setResidualToPCMat(pRes_, "p", PCMat);

    // ERes
    this->setResidualToPCMat(eRes_, "e", PCMat);

    // phiRes
    this->setResidualToPCMat(phiRes_, "phi", PCMat);
}

} // End namespace Foam

// ************************************************************************* //
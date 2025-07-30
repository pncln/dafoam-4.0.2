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
      setResidualClassMemberScalar(T, dimensionSet(1, -1, -3, 0, 0, 0, 0)),
      setResidualClassMemberPhi(phi),
      fvSource_(const_cast<volVectorField&>(
          mesh_.thisDb().lookupObject<volVectorField>("fvSource"))),
      fvSourceEnergy_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("fvSourceEnergy"))),
      thermo_(const_cast<psiThermo&>(
          mesh_.thisDb().lookupObject<psiThermo>("thermophysicalProperties"))),
      e_(thermo_.he()),
      rho_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("rho"))),
      alphat_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("alphat"))),
      psi_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("thermo:psi"))),
     // dpdt_(const_cast<volScalarField&>(
       //   mesh_.thisDb().lookupObject<volScalarField>("dpdt"))),
      K_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("K"))),
      daTurb_(const_cast<DATurbulenceModel&>(daModel.getDATurbulenceModel())),
      // create pimpleControl
      pimple_(const_cast<fvMesh&>(mesh))
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
    TRes_.clear();
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

        p_, T_, U_, phi_, etc: State variables in OpenFOAM
    
    Output:
        URes_, pRes_, TRes_, phiRes_: residual field variables
    */

    label isPC = options.getLabel("isPC");

    word divUScheme = "div(phi,U)";
    word divHEScheme = "div(phi,e)";

// 8314.4700665  gas constant in OpenFOAM
    // src/OpenFOAM/global/constants/thermodynamic/thermodynamicConstants.H
    scalar RR = Foam::constant::thermodynamic::RR;
  // R = RR/molWeight
    // Foam::specie::R() function in src/thermophysicalModels/specie/specie/specieI.H
    dimensionedScalar R(
        "R1",
        dimensionSet(0, 2, -2, -1, 0, 0, 0),
        RR / molWeight_);
 
   dimensionedScalar Cp(
        "Cp1",
        dimensionSet(0, 2, -2, -1, 0, 0, 0),
        Cp_);
	
    if (e_.name() == "e")
    {
        e_ = Cp * T_ - T_ * R;
        e_.oldTime() = Cp * T_.oldTime() - T_.oldTime() * R;
        e_.oldTime().oldTime() = Cp * T_.oldTime().oldTime() - T_.oldTime().oldTime() * R;
    }
    else
    {
        e_ = Cp * T_;
        e_.oldTime() = Cp * T_.oldTime();
        e_.oldTime().oldTime() = Cp * T_.oldTime().oldTime();
    }
    e_.correctBoundaryConditions();

    K_ = 0.5 * magSqr(U_);
    K_.oldTime() = 0.5 * magSqr(U_.oldTime());
    K_.oldTime().oldTime() = 0.5 * magSqr(U_.oldTime().oldTime());

   // dpdt_ = fvc::ddt(p_);

    // NOTE: alphat is updated in the correctNut function in DATurbulenceModel child classes
}

void Foam::DAResidualSonicFoam::updateIntermediateVariables()
{
    
}


void DAResidualSonicFoam::correctBoundaryConditions()
{
    /* 
    Description:
        Update the boundary condition for all the states in the selected solver
    */

    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    T_.correctBoundaryConditions();
}

void DAResidualSonicFoam::calcPCMatWithFvMatrix(Mat PCMat)
{
    /* 
    Description:
        Calculate the diagonal block of the preconditioner matrix dRdWTPC using the fvMatrix
    */

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    PetscScalar val;

    dictionary normStateDict = daOption_.getAllOptions().subDict("normalizeStates");
    wordList normResDict = daOption_.getOption<wordList>("normalizeResiduals");

    // ******** U Residuals **********

    fvVectorMatrix UEqn(
        fvm::ddt(rho_, U_)
        + fvm::div(phi_, U_, "div(pc)")
        + daTurb_.divDevReff(U_)
        - fvSource_);

    // set the val before relaxing UEqn!

    scalar UScaling = 1.0;
    if (normStateDict.found("U"))
    {
        UScaling = normStateDict.getScalar("U");
    }
    scalar UResScaling = 1.0;

    // set diag
    forAll(U_, cellI)
    {
        if (normResDict.found("URes"))
        {
            UResScaling = mesh_.V()[cellI];
        }
        for (label i = 0; i < 3; i++)
        {
            PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("U", cellI, i);
            PetscInt colI = rowI;
            scalarField D = UEqn.D();
            scalar val1 = D[cellI] * UScaling / UResScaling;
            assignValueCheckAD(val, val1);
            MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
        }
    }

    // set lower/owner
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("URes"))
        {
            UResScaling = mesh_.V()[neighbourCellI];
        }

        for (label i = 0; i < 3; i++)
        {
            PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("U", neighbourCellI, i);
            PetscInt colI = daIndex_.getGlobalAdjointStateIndex("U", ownerCellI, i);
            scalar val1 = UEqn.lower()[faceI] * UScaling / UResScaling;
            assignValueCheckAD(val, val1);
            MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
        }
    }

    // set upper/neighbour
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("URes"))
        {
            UResScaling = mesh_.V()[ownerCellI];
        }

        for (label i = 0; i < 3; i++)
        {
            PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("U", ownerCellI, i);
            PetscInt colI = daIndex_.getGlobalAdjointStateIndex("U", neighbourCellI, i);
            scalar val1 = UEqn.upper()[faceI] * UScaling / UResScaling;
            assignValueCheckAD(val, val1);
            MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
        }
    }

    UEqn.relax();

    // ******** e Residuals **********
    volScalarField alphaEff("alphaEff", thermo_.alphaEff(alphat_));

  fvScalarMatrix EEqn(
    fvm::ddt(rho_, e_)
  + fvm::div(phi_, e_, "div(pc)")
  + fvc::ddt(rho_, K_)
  + fvc::div(phi_, K_)
  + fvc::div( // e constant
        fvc::absolute(phi_ / fvc::interpolate(rho_), U_),
        p_,
        "div(phiv,p)")
  - fvm::laplacian(alphaEff, e_)
  - fvSourceEnergy_);


    scalar TScaling = 1.0;
    if (normStateDict.found("T"))
    {
        TScaling = normStateDict.getScalar("T");
    }
    scalar TResScaling = 1.0;

    // set diag
    forAll(T_, cellI)
    {
        if (normResDict.found("TRes"))
        {
            TResScaling = mesh_.V()[cellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("T", cellI);
        PetscInt colI = rowI;
        scalarField D = EEqn.D();
        scalar val1 = D[cellI] * TScaling / TResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // set lower/owner
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("TRes"))
        {
            TResScaling = mesh_.V()[neighbourCellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("T", neighbourCellI);
        PetscInt colI = daIndex_.getGlobalAdjointStateIndex("T", ownerCellI);
        scalar val1 = EEqn.lower()[faceI] * TScaling / TResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }

    // set upper/neighbour
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("TRes"))
        {
            TResScaling = mesh_.V()[ownerCellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("T", ownerCellI);
        PetscInt colI = daIndex_.getGlobalAdjointStateIndex("T", neighbourCellI);
        scalar val1 = EEqn.upper()[faceI] * TScaling / TResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }

    EEqn.relax();

    // ******** p Residuals **********
    volScalarField rAU(1.0 / UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_ * rAU));
    //***************** NOTE *******************
    // constrainHbyA has been used since OpenFOAM-v1606; however, it may degrade the accuracy of derivatives
    // because constraining variables will create discontinuity. Here we have a option to use the old
    // implementation in OpenFOAM-3.0+ and before (no constraint for HbyA)
    autoPtr<volVectorField> HbyAPtr = nullptr;
    label useConstrainHbyA = daOption_.getOption<label>("useConstrainHbyA");
    if (useConstrainHbyA)
    {
        HbyAPtr.reset(new volVectorField(constrainHbyA(rAU * UEqn.H(), U_, p_)));
    }
    else
    {
        HbyAPtr.reset(new volVectorField("HbyA", U_));
        HbyAPtr() = rAU * UEqn.H();
    }
    volVectorField& HbyA = HbyAPtr();

    surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho_) * fvc::flux(HbyA));

    // NOTE: we don't support transonic = true

    fvScalarMatrix pEqn(
        fvm::ddt(psi_, p_)
        + fvc::div(phiHbyA)
        - fvm::laplacian(rhorAUf, p_));

    scalar pScaling = 1.0;
    if (normStateDict.found("p"))
    {
        pScaling = normStateDict.getScalar("p");
    }
    scalar pResScaling = 1.0;
    // set diag
    forAll(p_, cellI)
    {
        if (normResDict.found("pRes"))
        {
            pResScaling = mesh_.V()[cellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p", cellI);
        PetscInt colI = rowI;
        scalarField D = pEqn.D();
        scalar val1 = D[cellI] * pScaling / pResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // set lower/owner
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("pRes"))
        {
            pResScaling = mesh_.V()[neighbourCellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p", neighbourCellI);
        PetscInt colI = daIndex_.getGlobalAdjointStateIndex("p", ownerCellI);
        scalar val1 = pEqn.lower()[faceI] * pScaling / pResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }

    // set upper/neighbour
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        if (normResDict.found("pRes"))
        {
            pResScaling = mesh_.V()[ownerCellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p", ownerCellI);
        PetscInt colI = daIndex_.getGlobalAdjointStateIndex("p", neighbourCellI);
        scalar val1 = pEqn.upper()[faceI] * pScaling / pResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }
}

} // End namespace Foam

// ************************************************************************* //

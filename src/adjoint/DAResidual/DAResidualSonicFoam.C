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
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

DAResidualSonicFoam::DAResidualSonicFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAResidual(modelType, mesh, daOption, daModel, daIndex),
      // initialize and register state variables and their residuals
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
      rho_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("rho"))),
      alphat_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("alphat"))),
      psi_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("thermo:psi"))),
      K_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("K"))),
      daTurb_(const_cast<DATurbulenceModel&>(daModel.getDATurbulenceModel())),
      pimple_(const_cast<fvMesh&>(mesh))
{
    // initialize fvSource
    const dictionary& allOptions = daOption.getAllOptions();
    if (allOptions.subDict("fvSource").toc().size() != 0)
    {
        hasFvSource_ = 1;
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
    */

    label isPC = options.getLabel("isPC");

    word divUScheme = "div(phi,U)";
    word divEScheme = "div(phi,e)";
    
    if (isPC)
    {
        divUScheme = "div(pc)";
        divEScheme = "div(pc)";
    }

    // ******** U Residuals **********
    if (hasFvSource_)
    {
        DAFvSource& daFvSource(const_cast<DAFvSource&>(
            mesh_.thisDb().lookupObject<DAFvSource>("DAFvSource")));
        daFvSource.calcFvSource(fvSource_);
        fvSourceEnergy_ = fvSource_ & U_;
    }

    fvVectorMatrix UEqn(
        fvm::ddt(rho_, U_)
        + fvm::div(phi_, U_, divUScheme)
        + daTurb_.divDevRhoReff(U_)
        - fvSource_);
    UEqn.relax(1.0);

    URes_ = (UEqn & U_) + fvc::grad(p_);
    normalizeResiduals(URes);

    // ******** T Residuals (computed from energy equation) **********
    // Get a reference to thermo's internal energy field
    volScalarField& e = thermo_.he();
    
    volScalarField alphaEff("alphaEff", thermo_.alphaEff(alphat_));

    fvScalarMatrix eEqn(
        fvm::ddt(rho_, e)
        + fvm::div(phi_, e, divEScheme)
        + fvc::ddt(rho_, K_)
        + fvc::div(phi_, K_)
        + fvc::div(
            fvc::absolute(phi_ / fvc::interpolate(rho_), U_),
            p_,
            "div(phiv,p)")
        - fvm::laplacian(alphaEff, e));

    if (hasFvSource_)
    {
        eEqn -= fvSourceEnergy_;
    }

    // The residual for T equation is computed from the energy equation
    TRes_ = eEqn & e;
    normalizeResiduals(TRes);

    // ******** p Residuals **********
    // Need to create a separate UEqn for pressure equation
    fvVectorMatrix UEqnP(
        fvm::ddt(rho_, U_)
        + fvm::div(phi_, U_)
        + daTurb_.divDevRhoReff(U_)
        - fvSource_);

    volScalarField rAU(1.0 / UEqnP.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_ * rAU));
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_ * rAU));
    
    volVectorField HbyA("HbyA", U_);
    label useConstrainHbyA = daOption_.getOption<label>("useConstrainHbyA");
    if (useConstrainHbyA)
    {
        HbyA = constrainHbyA(rAU * UEqnP.H(), U_, p_);
    }
    else
    {
        HbyA = rAU * UEqnP.H();
    }

    surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho_) * fvc::flux(HbyA));

    // Create phid for compressible pressure equation

    fvScalarMatrix pEqn(
        fvm::ddt(psi_, p_)
        + fvc::div(phiHbyA)
        - fvm::laplacian(rhorAUf, p_));

    pRes_ = pEqn & p_;
    normalizeResiduals(pRes);

    // ******** phi Residuals **********
    phiRes_ = phiHbyA + pEqn.flux() - phi_;
    normalizePhiResiduals(phiRes);
}

void DAResidualSonicFoam::updateIntermediateVariables()
{
    /*
    Description:
        Update intermediate variables that depend on state variables
    */
    
    // Sync thermo's T field with our T state
    thermo_.T() = T_;
    
    // Now correct thermo properties
    thermo_.correct();
    
    // Update density from equation of state
    rho_ = thermo_.rho();
    
    // Update compressibility
    // psi_ aliases thermo_.psi(); avoid self-assignment
    
    // Update kinetic energy
    K_ = 0.5 * magSqr(U_);
    
    // Update turbulence variables
    daTurb_.correctNut();
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
    
    // Update thermo after boundary corrections
    thermo_.T() = T_;
    thermo_.correct();
}

void DAResidualSonicFoam::calcPCMatWithFvMatrix(Mat PCMat)
{
    /* 
    Description:
        Calculate the diagonal block of the preconditioner matrix dRdWT using the fvMatrix.
        NOTE: Since we solve for T but the energy equation is for e, we need to
        properly account for the transformation de/dT in the preconditioner.
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
        + daTurb_.divDevRhoReff(U_)
        - fvSource_);

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

    // ******** T Residuals (from energy equation) **********
    volScalarField alphaEff("alphaEff", thermo_.alphaEff(alphat_));

    // Get reference to thermo's internal energy
    volScalarField& e = thermo_.he();

    fvScalarMatrix eEqn(
        fvm::ddt(rho_, e)
        + fvm::div(phi_, e, "div(pc)")
        + fvc::ddt(rho_, K_)
        + fvc::div(phi_, K_)
        + fvc::div(
            fvc::absolute(phi_ / fvc::interpolate(rho_), U_),
            p_,
            "div(phiv,p)")
        - fvm::laplacian(alphaEff, e));
    
    if (hasFvSource_)
    {
        eEqn -= fvSourceEnergy_;
    }

    // Get de/dT for the transformation (Cv for constant volume)
    volScalarField Cv = thermo_.Cv();
    
    scalar TScaling = 1.0;
    if (normStateDict.found("T"))
    {
        TScaling = normStateDict.getScalar("T");
    }
    scalar TResScaling = 1.0;

    // set diag - need to multiply by de/dT = Cv
    forAll(T_, cellI)
    {
        if (normResDict.found("TRes"))
        {
            TResScaling = mesh_.V()[cellI];
        }

        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("T", cellI);
        PetscInt colI = rowI;
        scalarField D = eEqn.D();
        scalar val1 = D[cellI] * Cv[cellI] * TScaling / TResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // set lower/owner - multiply by Cv at owner cell
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
        scalar val1 = eEqn.lower()[faceI] * Cv[ownerCellI] * TScaling / TResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }

    // set upper/neighbour - multiply by Cv at neighbour cell
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
        scalar val1 = eEqn.upper()[faceI] * Cv[neighbourCellI] * TScaling / TResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }

    eEqn.relax();

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
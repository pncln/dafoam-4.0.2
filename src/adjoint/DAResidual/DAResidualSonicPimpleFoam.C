/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAResidualSonicPimpleFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAResidualSonicPimpleFoam, 0);
addToRunTimeSelectionTable(DAResidual, DAResidualSonicPimpleFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAResidualSonicPimpleFoam::DAResidualSonicPimpleFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAResidual(modelType, mesh, daOption, daModel, daIndex),
      // initialize and register state variables and their residuals, we use macros defined in macroFunctions.H
      setResidualClassMemberVector(U, dimensionSet(1, -2, -2, 0, 0, 0, 0)),        // kg/(m2*s2)
      setResidualClassMemberScalar(p, dimensionSet(1, -1, -3, 0, 0, 0, 0)),        // kg/(m*s3)
      setResidualClassMemberScalar(rho, dimensionSet(1, -3, -1, 0, 0, 0, 0)),      // kg/(m3*s)
      setResidualClassMemberScalar(T, dimensionSet(0, 0, -1, 1, 0, 0, 0)),         // K/s
      setResidualClassMemberPhi(phi),
      fvSource_(const_cast<volVectorField&>(
          mesh_.thisDb().lookupObject<volVectorField>("fvSource"))),
      fvSourceEnergy_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("fvSourceEnergy"))),
      daTurb_(const_cast<DATurbulenceModel&>(daModel.getDATurbulenceModel())),
      thermo_(const_cast<fluidThermo&>(mesh_.thisDb().lookupObject<fluidThermo>("thermophysicalProperties"))),
      // create pimpleControl
      pimple_(const_cast<fvMesh&>(mesh)),
      pressureControl_(const_cast<pressureControl&>(mesh_.thisDb().lookupObject<pressureControl>("pressureControl"))),
      shockSensor_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("shockSensor")))
{
    // initialize fvSource
    const dictionary& allOptions = daOption.getAllOptions();
    if (allOptions.subDict("fvSource").toc().size() != 0)
    {
        hasFvSource_ = 1;
    }

    // initialize energy source term if needed
    if (allOptions.subDict("fvSource").found("energy"))
    {
        hasEnergySource_ = 1;
    }

    // initialize shock capturing options
    if (allOptions.found("shockCapturing"))
    {
        const dictionary& shockDict = allOptions.subDict("shockCapturing");
        shockCapturingScheme_ = shockDict.lookupOrDefault<word>("scheme", "artificialViscosity");
        artificialViscosityCoeff_ = shockDict.lookupOrDefault<scalar>("coefficient", 0.1);
    }
    else
    {
        shockCapturingScheme_ = "artificialViscosity";
        artificialViscosityCoeff_ = 0.1;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DAResidualSonicPimpleFoam::clear()
{
    /*
    Description:
        Clear all members to avoid memory leak because we will initalize 
        multiple objects of DAResidual. Here we need to delete all members
        in the parent and child classes
    */
    URes_.clear();
    pRes_.clear();
    rhoRes_.clear();
    TRes_.clear();
    phiRes_.clear();
}

void DAResidualSonicPimpleFoam::calcResiduals(const dictionary& options)
{
    /*
    Description:
        This is the function to compute residuals for compressible supersonic flow.
    
    Input:
        options.isPC: 1 means computing residuals for preconditioner matrix.
        This essentially use 1st order scheme for discretization. The default
        is 0 (use a higher order scheme)

        options.isRef: 1 means computing residuals for reference case. This is
        typically used when validating gradient accuracy with finite difference.

        options.getPC: whether to compute the preconditioning matrix dRdW.
        
        options.zeroResiduals: 1 means set all residuals to zero.

        options.updateState: whether to update state. This is useful when 
        computing partial derivatives with respect to the design variables.
        If options.updateState == 1, we update the state variables. If 
        options.updateState == 0, we don't update the state variables.
        However, we update intermediate variables that are related to state
        variables, e.g., k and omega variables for nut in turbulence models.
    */

    word divUScheme = "div(phi,U)";
    word divTScheme = "div(phi,T)";
    word divPhiScheme = "div(phi)";

    label isPC = 0;

    if (!options.readIfPresent<label>("isPC", isPC))
    {
        isPC = 0;
    }

    if (isPC)
    {
        divUScheme = "div(pc)";
        divTScheme = "div(pc)";
        divPhiScheme = "div(pc)";
    }

    // We dont support zeroResiduals and updateState options for compressible flow yet
    label zeroResiduals = 0;
    options.readIfPresent<label>("zeroResiduals", zeroResiduals);

    label updateState = 1;
    options.readIfPresent<label>("updateState", updateState);

    // if zeroResiduals == 1, we set all the residuals to zero
    if (zeroResiduals)
    {
        URes_ = dimensionedVector("URes", URes_.dimensions(), vector::zero);
        pRes_ = dimensionedScalar("pRes", pRes_.dimensions(), 0.0);
        rhoRes_ = dimensionedScalar("rhoRes", rhoRes_.dimensions(), 0.0);
        TRes_ = dimensionedScalar("TRes", TRes_.dimensions(), 0.0);
        phiRes_ = dimensionedScalar("phiRes", phiRes_.dimensions(), 0.0);
        return;
    }

    // calculate shock sensor for artificial viscosity
    this->calcShockSensor();

    // bound thermodynamic properties
    this->boundThermodynamicProperties();

    // Continuity (Density) equation
    fvScalarMatrix rhoEqn(
        fvm::ddt(rho_)
      + fvc::div(phi_)
    );

    rhoRes_ = rhoEqn & rho_;

    // Momentum equation with compressible effects and shock capturing
    volScalarField muEff("muEff", daTurb_.muEff());
    
    // Add artificial viscosity for shock capturing
    if (shockCapturingScheme_ == "artificialViscosity")
    {
        volScalarField cellVolume = mesh_.V();
        volScalarField cellSize = pow(cellVolume, 1.0/3.0);
        volScalarField artificialViscosity = 
            artificialViscosityCoeff_ * shockSensor_ * cellSize * rho_ * mag(U_);
        muEff += artificialViscosity;
    }

    fvVectorMatrix UEqn(
        fvm::ddt(rho_, U_)
      + fvm::div(phi_, U_, divUScheme)
      + daTurb_.divDevRhoReff(U_)
      - fvSource_
    );

    URes_ = (UEqn & U_) + fvc::grad(p_);

    // Energy equation including kinetic energy effects
    volScalarField& he = thermo_.he();
    volScalarField alphaEff("alphaEff", daTurb_.alphaEff());
    
    // Add artificial thermal diffusivity for shock capturing
    if (shockCapturingScheme_ == "artificialViscosity")
    {
        volScalarField cellVolume = mesh_.V();
        volScalarField cellSize = pow(cellVolume, 1.0/3.0);
        volScalarField artificialAlpha = 
            artificialViscosityCoeff_ * shockSensor_ * cellSize * rho_ * mag(U_) * thermo_.Cp();
        alphaEff += artificialAlpha;
    }

    fvScalarMatrix EEqn(
        fvm::ddt(rho_, he)
      + fvm::div(phi_, he, divTScheme)
      + fvc::ddt(rho_, 0.5*magSqr(U_))
      + fvc::div(phi_, 0.5*magSqr(U_))
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi_/fvc::interpolate(rho_), U_),
                p_,
                "div(phiv,p)"
            )
          : -fvc::ddt(p_)
        )
      - fvm::laplacian(alphaEff, he)
      - fvSourceEnergy_
    );

    TRes_ = EEqn & he;

    // Pressure equation with acoustic effects
    tmp<volScalarField> trAU(1.0/UEqn.A());
    const volScalarField& rAU = trAU();

    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, p_));
    surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho_)*fvc::flux(HbyA));

    // Get compressibility
    tmp<volScalarField> tpsi(thermo_.psi());
    const volScalarField& psi = tpsi();

    fvScalarMatrix pEqn(
        fvm::ddt(psi, p_)
      + fvm::div(phiHbyA, divPhiScheme)
      - fvm::laplacian(rho_*rAU, p_)
    );

    pRes_ = pEqn & p_;

    // Face flux residual
    phiRes_ = phiHbyA + pEqn.flux() - phi_;
}

void DAResidualSonicPimpleFoam::updateIntermediateVariables()
{
    /*
    Description:
        Update the intermediate variables that depend on the state variables
    */

    // update density from equation of state
    rho_ = thermo_.rho();

    // update boundary conditions for rho
    rho_.correctBoundaryConditions();

    // update thermodynamic properties
    thermo_.correct();

    // update turbulence variables
    daTurb_.correct();
}

void DAResidualSonicPimpleFoam::correctBoundaryConditions()
{
    /*
    Description:
        Update the boundary condition for all the states in the selected solver
    */

    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    T_.correctBoundaryConditions();
    phi_.correctBoundaryConditions();
}

void DAResidualSonicPimpleFoam::calcPCMatWithFvMatrix(Mat PCMat)
{
    /*
    Description:
        Calculate the preconditioner matrix dRdW for compressible flow.
        This is a simplified diagonal preconditioning approach.
    */

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();

    PetscScalar val = 0.0;
    
    // Get scaling factors from normalization
    const dictionary& allOptions = daOption_.getAllOptions();
    const dictionary& normStateDict = allOptions.subDict("normalizeStates");
    const dictionary& normResDict = allOptions.subDict("normalizeResiduals");
    
    scalar UScaling = 1.0;
    scalar pScaling = 1.0;
    scalar rhoScaling = 1.0;
    scalar TScaling = 1.0;
    scalar phiScaling = 1.0;
    
    scalar UResScaling = 1.0;
    scalar pResScaling = 1.0;
    scalar rhoResScaling = 1.0;
    scalar TResScaling = 1.0;
    scalar phiResScaling = 1.0;

    if (normStateDict.found("U"))
    {
        UScaling = normStateDict.getScalar("U");
    }
    if (normStateDict.found("p"))
    {
        pScaling = normStateDict.getScalar("p");
    }
    if (normStateDict.found("rho"))
    {
        rhoScaling = normStateDict.getScalar("rho");
    }
    if (normStateDict.found("T"))
    {
        TScaling = normStateDict.getScalar("T");
    }
    if (normStateDict.found("phi"))
    {
        phiScaling = normStateDict.getScalar("phi");
    }

    // Calculate simplified diagonal preconditioner entries
    // This is a basic implementation - could be enhanced with coupling terms

    // U diagonal entries
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
            scalar val1 = 1.0 * UScaling / UResScaling;  // Simplified diagonal
            assignValueCheckAD(val, val1);
            MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
        }
    }

    // p diagonal entries
    forAll(p_, cellI)
    {
        if (normResDict.found("pRes"))
        {
            pResScaling = mesh_.V()[cellI];
        }
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p", cellI);
        PetscInt colI = rowI;
        scalar val1 = 1.0 * pScaling / pResScaling;  // Simplified diagonal
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // rho diagonal entries
    forAll(rho_, cellI)
    {
        if (normResDict.found("rhoRes"))
        {
            rhoResScaling = mesh_.V()[cellI];
        }
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("rho", cellI);
        PetscInt colI = rowI;
        scalar val1 = 1.0 * rhoScaling / rhoResScaling;  // Simplified diagonal
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // T diagonal entries
    forAll(T_, cellI)
    {
        if (normResDict.found("TRes"))
        {
            TResScaling = mesh_.V()[cellI];
        }
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("T", cellI);
        PetscInt colI = rowI;
        scalar val1 = 1.0 * TScaling / TResScaling;  // Simplified diagonal
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // phi diagonal entries
    forAll(phi_, faceI)
    {
        if (normResDict.found("phiRes"))
        {
            phiResScaling = 1.0;
        }
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("phi", faceI);
        PetscInt colI = rowI;
        scalar val1 = 1.0 * phiScaling / phiResScaling;  // Simplified diagonal
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }
}

void DAResidualSonicPimpleFoam::calcShockSensor()
{
    /*
    Description:
        Calculate shock sensor based on pressure gradient magnitude
    */
    
    // Calculate pressure gradient magnitude
    volVectorField gradP = fvc::grad(p_);
    volScalarField magGradP = mag(gradP);
    
    // Normalize by local pressure to get relative gradient
    volScalarField relativeGradP = magGradP / (p_ + dimensionedScalar("small", p_.dimensions(), SMALL));
    
    // Apply simple averaging to smooth oscillations (instead of fvc::smooth)
    shockSensor_ = relativeGradP;
    
    // Apply a simple 3x3 cell averaging for smoothing
    forAll(shockSensor_, cellI)
    {
        scalar avgSensor = relativeGradP[cellI];
        label nNeighbors = 1;
        
        const labelList& cellCells = mesh_.cellCells()[cellI];
        forAll(cellCells, neighborI)
        {
            avgSensor += relativeGradP[cellCells[neighborI]];
            nNeighbors++;
        }
        
        shockSensor_[cellI] = avgSensor / nNeighbors;
    }
    
    // Update boundary conditions
    shockSensor_.correctBoundaryConditions();
    
    // Limit shock sensor to reasonable values
    shockSensor_ = min(shockSensor_, dimensionedScalar("maxSensor", dimless, 10.0));
}

void DAResidualSonicPimpleFoam::applyArtificialViscosity()
{
    /*
    Description:
        Apply artificial viscosity for shock capturing
        Note: This function is called during residual calculation
    */
    
    // The artificial viscosity is applied directly in calcResiduals()
    // This function is kept for potential future enhancements
}

void DAResidualSonicPimpleFoam::boundThermodynamicProperties()
{
    /*
    Description:
        Bound density, pressure and temperature to physical values
    */
    
    // Bound density to positive values
    rho_ = max(rho_, dimensionedScalar("rhoMin", rho_.dimensions(), 1e-8));
    rho_ = min(rho_, dimensionedScalar("rhoMax", rho_.dimensions(), 1000.0));
    
    // Bound pressure to positive values (use a minimum pressure value)
    dimensionedScalar pMin("pMin", p_.dimensions(), 1e-8);
    p_ = max(p_, pMin);
    
    // Bound temperature
    T_ = max(T_, dimensionedScalar("TMin", T_.dimensions(), 50.0));
    T_ = min(T_, dimensionedScalar("TMax", T_.dimensions(), 5000.0));
}

} // End namespace Foam

// ************************************************************************* //
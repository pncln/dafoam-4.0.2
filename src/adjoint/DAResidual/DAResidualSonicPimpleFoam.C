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
    // For compressible flow: dynamic viscosity = kinematic viscosity * density
    volScalarField muEff("muEff", daTurb_.nuEff() * rho_);
    
    // Add artificial viscosity for shock capturing
    if (shockCapturingScheme_ == "artificialViscosity")
    {
        // Create proper volScalarField for cell volume
        volScalarField cellVolume(
            IOobject(
                "cellVolume",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar("zero", dimVolume, 0.0));
        
        forAll(cellVolume, cellI)
        {
            cellVolume[cellI] = mesh_.V()[cellI];
        }
        
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
        volScalarField cellVolume(
            IOobject(
                "cellVolume2",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar("zero", dimVolume, 0.0));
        
        forAll(cellVolume, cellI)
        {
            cellVolume[cellI] = mesh_.V()[cellI];
        }
        
        volScalarField cellSize = pow(cellVolume, 1.0/3.0);
        volScalarField artificialAlpha = 
            artificialViscosityCoeff_ * shockSensor_ * cellSize * rho_ * mag(U_) * thermo_.Cp();
        alphaEff += artificialAlpha;
    }

    // Calculate kinetic energy terms properly
    volScalarField K = 0.5 * magSqr(U_);
    volScalarField rhoK = rho_ * K;
    
    fvScalarMatrix EEqn(
        fvm::ddt(rho_, he)
      + fvm::div(phi_, he, divTScheme)
      + fvc::ddt(rhoK)
      + fvc::div(phi_, K)
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
      + fvc::div(phiHbyA)
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
    daTurb_.correct(0);
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
        Calculate the preconditioner matrix dRdW for compressible supersonic flow.
        This includes proper momentum, pressure, density, and energy coupling.
    */

    const labelList& owner = mesh_.owner();
    const labelList& neighbour = mesh_.neighbour();
    const dictionary& allOptions = daOption_.getAllOptions();
    
    PetscScalar val = 0.0;
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

    // Get normalization factors if they exist
    if (allOptions.found("normalizeStates"))
    {
        const dictionary& normStateDict = allOptions.subDict("normalizeStates");
        if (normStateDict.found("U")) UScaling = normStateDict.getScalar("U");
        if (normStateDict.found("p")) pScaling = normStateDict.getScalar("p");
        if (normStateDict.found("rho")) rhoScaling = normStateDict.getScalar("rho");
        if (normStateDict.found("T")) TScaling = normStateDict.getScalar("T");
        if (normStateDict.found("phi")) phiScaling = normStateDict.getScalar("phi");
    }

    if (allOptions.found("normalizeResiduals"))
    {
        const dictionary& normResDict = allOptions.subDict("normalizeResiduals");
        // Residual scaling will be handled per cell if needed
    }

    // ******** Momentum equation matrix (URes) **********
    volScalarField muEff("muEff", daTurb_.nuEff() * rho_);
    
    // Add artificial viscosity for shock capturing
    if (shockCapturingScheme_ == "artificialViscosity")
    {
        volScalarField cellVolume(
            IOobject("cellVolume", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_, dimensionedScalar("zero", dimVolume, 0.0));
        forAll(cellVolume, cellI) { cellVolume[cellI] = mesh_.V()[cellI]; }
        
        volScalarField cellSize = pow(cellVolume, 1.0/3.0);
        muEff += artificialViscosityCoeff_ * shockSensor_ * cellSize * rho_ * mag(U_);
    }

    fvVectorMatrix UEqn(
        fvm::ddt(rho_, U_)
      + fvm::div(phi_, U_, "div(pc)")
      + daTurb_.divDevRhoReff(U_)
      - fvSource_);
    UEqn.relax(1.0);

    // Set diagonal entries for momentum
    forAll(U_, cellI)
    {
        UResScaling = mesh_.V()[cellI];
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

    // Set off-diagonal entries for momentum
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        UResScaling = mesh_.V()[neighbourCellI];

        for (label i = 0; i < 3; i++)
        {
            PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("U", neighbourCellI, i);
            PetscInt colI = daIndex_.getGlobalAdjointStateIndex("U", ownerCellI, i);
            scalar val1 = UEqn.lower()[faceI] * UScaling / UResScaling;
            assignValueCheckAD(val, val1);
            MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
        }
    }

    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        UResScaling = mesh_.V()[ownerCellI];

        for (label i = 0; i < 3; i++)
        {
            PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("U", ownerCellI, i);
            PetscInt colI = daIndex_.getGlobalAdjointStateIndex("U", neighbourCellI, i);
            scalar val1 = UEqn.upper()[faceI] * UScaling / UResScaling;
            assignValueCheckAD(val, val1);
            MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
        }
    }

    // ******** Pressure equation matrix (pRes) **********
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, p_));
    surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho_)*fvc::flux(HbyA));
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_*rAU));

    fvScalarMatrix pEqn(
        fvm::ddt(thermo_.psi(), p_)
      + fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, p_));

    // Set diagonal entries for pressure
    forAll(p_, cellI)
    {
        pResScaling = mesh_.V()[cellI];
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p", cellI);
        PetscInt colI = rowI;
        scalarField D = pEqn.D();
        scalar val1 = D[cellI] * pScaling / pResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // Set off-diagonal entries for pressure
    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        pResScaling = mesh_.V()[neighbourCellI];
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p", neighbourCellI);
        PetscInt colI = daIndex_.getGlobalAdjointStateIndex("p", ownerCellI);
        scalar val1 = pEqn.lower()[faceI] * pScaling / pResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }

    for (label faceI = 0; faceI < daIndex_.nLocalInternalFaces; faceI++)
    {
        label ownerCellI = owner[faceI];
        label neighbourCellI = neighbour[faceI];

        pResScaling = mesh_.V()[ownerCellI];
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("p", ownerCellI);
        PetscInt colI = daIndex_.getGlobalAdjointStateIndex("p", neighbourCellI);
        scalar val1 = pEqn.upper()[faceI] * pScaling / pResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &colI, 1, &rowI, &val, INSERT_VALUES);
    }

    // ******** Energy/Temperature equation matrix (TRes) **********
    volScalarField& he = thermo_.he();
    volScalarField alphaEff("alphaEff", daTurb_.alphaEff());
    
    // Add artificial thermal diffusivity
    if (shockCapturingScheme_ == "artificialViscosity")
    {
        volScalarField cellVolume(
            IOobject("cellVolume2", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_, dimensionedScalar("zero", dimVolume, 0.0));
        forAll(cellVolume, cellI) { cellVolume[cellI] = mesh_.V()[cellI]; }
        
        volScalarField cellSize = pow(cellVolume, 1.0/3.0);
        alphaEff += artificialViscosityCoeff_ * shockSensor_ * cellSize * rho_ * mag(U_) * thermo_.Cp();
    }

    volScalarField K = 0.5 * magSqr(U_);
    volScalarField rhoK = rho_ * K;
    
    fvScalarMatrix EEqn(
        fvm::ddt(rho_, he)
      + fvm::div(phi_, he, "div(pc)")
      + fvc::ddt(rhoK)
      + fvc::div(phi_, K)
      - fvm::laplacian(alphaEff, he)
      - fvSourceEnergy_);

    // Set diagonal entries for energy/temperature
    forAll(T_, cellI)
    {
        TResScaling = mesh_.V()[cellI];
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("T", cellI);
        PetscInt colI = rowI;
        scalarField D = EEqn.D();
        scalar val1 = D[cellI] * TScaling / TResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // ******** Density equation matrix (rhoRes) **********
    fvScalarMatrix rhoEqn(fvm::ddt(rho_) + fvc::div(phi_));

    forAll(rho_, cellI)
    {
        rhoResScaling = mesh_.V()[cellI];
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("rho", cellI);
        PetscInt colI = rowI;
        scalarField D = rhoEqn.D();
        scalar val1 = D[cellI] * rhoScaling / rhoResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }

    // ******** Flux equation (phiRes) - diagonal **********
    forAll(phi_, faceI)
    {
        phiResScaling = 1.0;
        PetscInt rowI = daIndex_.getGlobalAdjointStateIndex("phi", faceI);
        PetscInt colI = rowI;
        scalar val1 = 1.0 * phiScaling / phiResScaling;
        assignValueCheckAD(val, val1);
        MatSetValues(PCMat, 1, &rowI, 1, &colI, &val, INSERT_VALUES);
    }
}

void DAResidualSonicPimpleFoam::calcShockSensor()
{
    /*
    Description:
        Calculate advanced shock sensor for supersonic flows using multiple indicators:
        1. Pressure gradient magnitude (primary indicator)
        2. Velocity divergence (expansion/compression detection)
        3. Density gradient (density jumps across shocks)
        4. Mach number gradient (transonic/supersonic transitions)
    */
    
    // Primary shock indicator: normalized pressure gradient
    volVectorField gradP = fvc::grad(p_);
    volScalarField magGradP = mag(gradP);
    volScalarField pressureIndicator = magGradP / 
        (p_ + dimensionedScalar("small", p_.dimensions(), SMALL));
    
    // Secondary indicator: velocity divergence (detects expansion fans and compression)
    volScalarField divU = fvc::div(U_);
    volScalarField compressionIndicator = mag(divU) * mag(U_) / 
        (mag(U_) + dimensionedScalar("smallU", U_.dimensions(), SMALL));
    
    // Tertiary indicator: density gradient (captures density jumps)
    volVectorField gradRho = fvc::grad(rho_);
    volScalarField magGradRho = mag(gradRho);
    volScalarField densityIndicator = magGradRho / 
        (rho_ + dimensionedScalar("smallRho", rho_.dimensions(), SMALL));
    
    // Quaternary indicator: Mach number gradient for transonic regions
    volScalarField c = sqrt(thermo_.Cp()/thermo_.Cv() * p_ / rho_); // Speed of sound
    volScalarField Mach = mag(U_) / c;
    volVectorField gradMach = fvc::grad(Mach);
    volScalarField machIndicator = mag(gradMach);
    
    // Combine indicators with weighted average
    // Higher weights for primary indicators in supersonic flows
    scalar wP = 0.5;    // Pressure gradient weight
    scalar wDiv = 0.2;  // Velocity divergence weight  
    scalar wRho = 0.2;  // Density gradient weight
    scalar wMach = 0.1; // Mach gradient weight
    
    volScalarField combinedSensor = 
        wP * pressureIndicator + 
        wDiv * compressionIndicator + 
        wRho * densityIndicator + 
        wMach * machIndicator;
    
    // Apply advanced multi-pass smoothing to reduce numerical noise
    // Pass 1: Weighted average with immediate neighbors
    volScalarField smoothedSensor = combinedSensor;
    forAll(smoothedSensor, cellI)
    {
        scalar weightedSum = combinedSensor[cellI];
        scalar totalWeight = 1.0;
        
        const labelList& cellCells = mesh_.cellCells()[cellI];
        forAll(cellCells, neighborI)
        {
            label neighCellI = cellCells[neighborI];
            scalar weight = 0.5; // Neighbor weight
            weightedSum += weight * combinedSensor[neighCellI];
            totalWeight += weight;
        }
        
        smoothedSensor[cellI] = weightedSum / totalWeight;
    }
    
    // Pass 2: Apply TVD limiter to prevent spurious oscillations
    forAll(smoothedSensor, cellI)
    {
        scalar sensorValue = smoothedSensor[cellI];
        scalar maxNeighbor = sensorValue;
        scalar minNeighbor = sensorValue;
        
        const labelList& cellCells = mesh_.cellCells()[cellI];
        forAll(cellCells, neighborI)
        {
            label neighCellI = cellCells[neighborI];
            maxNeighbor = max(maxNeighbor, smoothedSensor[neighCellI]);
            minNeighbor = min(minNeighbor, smoothedSensor[neighCellI]);
        }
        
        // TVD limiting: constrain sensor between min and max of neighbors
        scalar limitedValue = max(minNeighbor, min(maxNeighbor, sensorValue));
        smoothedSensor[cellI] = limitedValue;
    }
    
    // Update boundary conditions
    smoothedSensor.correctBoundaryConditions();
    
    // Apply nonlinear scaling to enhance shock detection
    // Use tanh function for smooth activation near shocks
    scalar threshold = 0.1; // Threshold for shock detection
    scalar steepness = 10.0; // Steepness of activation function
    
    forAll(shockSensor_, cellI)
    {
        scalar rawSensor = smoothedSensor[cellI];
        // Smooth step function: 0.5 * (1 + tanh(steepness * (sensor - threshold)))
        scalar activatedSensor = 0.5 * (1.0 + tanh(steepness * (rawSensor - threshold)));
        shockSensor_[cellI] = activatedSensor;
    }
    
    // Final bounds to ensure physical values
    shockSensor_ = max(shockSensor_, dimensionedScalar("minSensor", dimless, 0.0));
    shockSensor_ = min(shockSensor_, dimensionedScalar("maxSensor", dimless, 1.0));
    
    // Update boundary conditions one final time
    shockSensor_.correctBoundaryConditions();
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
Info << "Reading thermophysical properties\n" << endl;

pThermoPtr_.reset(fluidThermo::New(mesh).ptr());
fluidThermo& pThermo = pThermoPtr_();
pThermo.validate(args.executable(), "h", "e");

volScalarField& p = pThermo.p();
pPtr_.reset(&p);

const volScalarField& T = pThermo.T();
const volScalarField& psi = pThermo.psi();

Info << "Reading field U\n" << endl;
UPtr_.reset(
    new volVectorField(
        IOobject(
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh));
volVectorField& U = UPtr_();

Info << "Creating density field rho\n" << endl;
rhoPtr_.reset(
    new volScalarField(
        IOobject(
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        pThermo.rho()));
volScalarField& rho = rhoPtr_();

Info << "Creating face flux field phi\n" << endl;
phiPtr_.reset(
    new surfaceScalarField(
        IOobject(
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        linearInterpolate(rho*U) & mesh.Sf()));
surfaceScalarField& phi = phiPtr_();

// Create pressure control for closed domains
pressureControlPtr_.reset(
    new pressureControl(p, rho, pimple.dict(), false));
pressureControl& pressureControl = pressureControlPtr_();

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model\n" << endl;
turbulencePtr_.reset(
    compressible::turbulenceModel::New(
        rho,
        U,
        phi,
        pThermo
    ).ptr());
compressible::turbulenceModel& turbulence = turbulencePtr_();

// Create source terms
Info << "Creating momentum source term\n" << endl;
fvSourcePtr_.reset(
    new volVectorField(
        IOobject(
            "fvSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("fvSource", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero),
        zeroGradientFvPatchField<vector>::typeName));
volVectorField& fvSource = fvSourcePtr_();

Info << "Creating energy source term\n" << endl;
fvSourceEnergyPtr_.reset(
    new volScalarField(
        IOobject(
            "fvSourceEnergy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("fvSourceEnergy", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchField<scalar>::typeName));
volScalarField& fvSourceEnergy = fvSourceEnergyPtr_();

#include "compressibleCreatePhi.H"
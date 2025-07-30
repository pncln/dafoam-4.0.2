Info<< "Reading thermophysical properties\n" << endl;

pThermoPtr_.reset(psiThermo::New(mesh).ptr());
psiThermo& thermo = pThermoPtr_();
thermo.validate(args.executable(), "e");

Info << "Reading field p\n" << endl;
pPtr_.reset(new volScalarField(
    IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh));
volScalarField& p = pPtr_();

rhoPtr_.reset(new volScalarField(
    IOobject("rho", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    thermo.rho()));
volScalarField& rho = rhoPtr_();

Info << "Reading field U\n" << endl;
UPtr_.reset(new volVectorField(
    IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh));
volVectorField& U = UPtr_();

Info<< "Creating field internal energy e\n" << endl;
ePtr_.reset(new volScalarField(
    IOobject("e", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("e", dimEnergy/dimMass, 0.0)));

// Create phi from compressibleCreatePhi.H
#include "compressibleCreatePhi.H"
// surfaceScalarField phi is created here

// actuator related stuff
Info << "Creating source term. " << endl;
fvSourcePtr_.reset(new volVectorField(
    IOobject("fvSource", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("fvSource", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero),
    zeroGradientFvPatchField<vector>::typeName));

fvSourceEnergyPtr_.reset(new volScalarField(
    IOobject("fvSourceEnergy", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("fvSourceEnergy", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchField<scalar>::typeName));

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
turbulencePtr_.reset(compressible::turbulenceModel::New(rho, U, phi, thermo));

Info << "Creating field kinetic energy K\n" << endl;
KPtr_.reset(new volScalarField("K", 0.5 * magSqr(U)));

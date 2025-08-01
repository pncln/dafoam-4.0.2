Info<< "Reading thermophysical properties\n" << endl;

pThermoPtr_.reset(psiThermo::New(mesh).ptr());
psiThermo& thermo = pThermoPtr_();
thermo.validate(args.executable(), "e");  // Changed from "h" to "e"

Info << "Reading field p\n" << endl;
pPtr_.reset(new volScalarField(
    IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh));
volScalarField& p = pPtr_();

Info << "Reading field T\n" << endl;
TPtr_.reset(new volScalarField(
    IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh));
volScalarField& T = TPtr_();
(void)T;  // Suppress unused variable warning

Info << "Reading field U\n" << endl;
UPtr_.reset(new volVectorField(
    IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh));
volVectorField& U = UPtr_();

// Initialize density from equation of state
rhoPtr_.reset(
    new volScalarField(
        IOobject(
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        thermo.rho()));
volScalarField& rho = rhoPtr_();

// Create internal energy field (not read from disk, calculated from T)
Info<< "Creating internal energy field e\n" << endl;
ePtr_.reset(
    new volScalarField(
        IOobject(
            "e",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        thermo.he()  // Initialize from thermo model
    ));
volScalarField& e = ePtr_();
(void)e;  // Suppress unused variable warning

#include "compressibleCreatePhiPython.H"
phi.write();

// Actuator disk source terms
Info << "Creating source terms\n" << endl;
fvSourcePtr_.reset(new volVectorField(
    IOobject("fvSource", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("fvSource", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)));

fvSourceEnergyPtr_.reset(new volScalarField(
    IOobject("fvSourceEnergy", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("fvSourceEnergy", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)));

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
turbulencePtr_.reset(compressible::turbulenceModel::New(rho, U, phi, thermo));

Info << "Creating field dpdt\n"
     << endl;

dpdtPtr_.reset(
    new volScalarField(
        IOobject(
            "dpdt",
            runTime.timeName(),
            mesh),
        mesh,
        dimensionedScalar("dpdt", p.dimensions() / dimTime, 0.0)));

Info << "Creating field kinetic energy K\n" << endl;
KPtr_.reset(new volScalarField("K", 0.5 * magSqr(U)));

// Initialize thermo fields
thermo.correct();
rho = thermo.rho();
rho.oldTime();

#include "createFvOptions.H"

Info<< "Reading thermophysical properties\n" << endl;

pThermoPtr_.reset(psiThermo::New(mesh).ptr());
psiThermo& thermo = pThermoPtr_();
thermo.validate(args.executable(), "e");

// Read pressure field separately (not from thermo)
Info << "Reading field p\n" << endl;
pPtr_.reset(
    new volScalarField(
        IOobject(
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh));
volScalarField& p = pPtr_();

// Read temperature field separately
Info << "Reading field T\n" << endl;
TPtr_.reset(
    new volScalarField(
        IOobject(
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh));
volScalarField& T = TPtr_();

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

// Create rho field  
rhoPtr_.reset(
    new volScalarField(
        IOobject(
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("rho", dimDensity, 0.0)));
volScalarField& rho = rhoPtr_();

// Make sure thermo's p and T fields match our fields
thermo.p() = p;
thermo.T() = T;

// Now update thermo properties
thermo.correct();

// Update rho from thermo
rho = thermo.rho();

#include "compressibleCreatePhiPython.H"
phi.write();

// Actuator disk source terms
Info << "Creating source terms\n" << endl;
fvSourcePtr_.reset(
    new volVectorField(
        IOobject(
            "fvSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("fvSource", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)));

fvSourceEnergyPtr_.reset(
    new volScalarField(
        IOobject(
            "fvSourceEnergy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("fvSourceEnergy", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)));

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
turbulencePtr_.reset(
    compressible::turbulenceModel::New(
        rho,
        U,
        phi,
        thermo));

Info << "Creating field dpdt\n" << endl;
dpdtPtr_.reset(
    new volScalarField(
        IOobject(
            "dpdt",
            runTime.timeName(),
            mesh),
        mesh,
        dimensionedScalar("dpdt", p.dimensions() / dimTime, 0.0)));

Info << "Creating field kinetic energy K\n" << endl;
KPtr_.reset(
    new volScalarField(
        "K",
        0.5 * magSqr(U)));

// Store old time
rho.oldTime();

#include "createFvOptions.H"

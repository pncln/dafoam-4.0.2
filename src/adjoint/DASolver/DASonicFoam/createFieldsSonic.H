Info<< "Reading thermophysical properties\n" << endl;

pThermoPtr_.reset(psiThermo::New(mesh).ptr());
psiThermo& thermo = pThermoPtr_();
thermo.validate(args.executable(), "e");  // Changed from "h" to "e"

volScalarField& p = thermo.p();

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

rhoPtr_.reset(
    new volScalarField(
        IOobject(
            "rho",
            runTime.timeName(),
            mesh),
        thermo.rho()));
volScalarField& rho = rhoPtr_();

volScalarField& e = thermo.he();

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
rho.oldTime();

#include "createFvOptions.H"

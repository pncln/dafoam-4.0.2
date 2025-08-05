Time& runTime = runTimePtr_();
// reset time to 0
runTime.setTime(0.0, 0);
fvMesh& mesh = meshPtr_();
pimpleControlDF& pimple = pimplePtr_();
psiThermo& thermo = pThermoPtr_();
volScalarField& p = *pPtr_;  // dereference the raw pointer
volScalarField& T = TPtr_();
volScalarField& rho = rhoPtr_();
volVectorField& U = UPtr_();
volScalarField& dpdt = dpdtPtr_();
surfaceScalarField& phi = phiPtr_();
volScalarField& K = KPtr_();
// Access thermo.he() directly without storing a reference
compressible::turbulenceModel& turbulence = turbulencePtr_();
scalar& cumulativeContErr = cumulativeContErr_;
const volScalarField& psi = thermo.psi();
const dictionary& allOptions = daOptionPtr_->getAllOptions();
volVectorField& fvSource = fvSourcePtr_();
volScalarField& fvSourceEnergy = fvSourceEnergyPtr_();
fv::options fvOptions(mesh);

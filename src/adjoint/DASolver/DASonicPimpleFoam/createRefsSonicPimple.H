fluidThermo& pThermo = pThermoPtr_();
volScalarField& p = pPtr_();
volVectorField& U = UPtr_();
volScalarField& rho = rhoPtr_();
surfaceScalarField& phi = phiPtr_();
compressible::turbulenceModel& turbulence = turbulencePtr_();
volVectorField& fvSource = fvSourcePtr_();
volScalarField& fvSourceEnergy = fvSourceEnergyPtr_();
pressureControl& pressureControl = pressureControlPtr_();
volScalarField& shockSensor = shockSensorPtr_();

const volScalarField& T = pThermo.T();
const volScalarField& psi = pThermo.psi();
const dimensionedScalar& pMin = pThermo.pMin();
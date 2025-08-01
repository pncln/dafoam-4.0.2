/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    sonicFoam

Group
    grpCompressibleSolvers

Description
    Transient solver for trans-sonic/supersonic, turbulent flow of a
    compressible gas.

\*---------------------------------------------------------------------------*/

#include "DASonicFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

defineTypeNameAndDebug(DASonicFoam, 0);
addToRunTimeSelectionTable(DASolver, DASonicFoam, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DASonicFoam::DASonicFoam(
    char* argsAll,
    PyObject* pyOptions)
    : DASolver(argsAll, pyOptions),
      pimplePtr_(nullptr),
      pThermoPtr_(nullptr),
      pPtr_(nullptr),				// check thisssss!!!!!!!!!!!!!! //
      rhoPtr_(nullptr),
      UPtr_(nullptr),
      phiPtr_(nullptr),
      KPtr_(nullptr),
      turbulencePtr_(nullptr),
      daTurbulenceModelPtr_(nullptr),
      daFvSourcePtr_(nullptr),
      fvSourcePtr_(nullptr),
      fvSourceEnergyPtr_(nullptr)
{
    // Create a dummy phi field early to avoid lookup errors
    // This will be replaced with the proper phi in initSolver()
    Time& runTime = runTimePtr_();
    fvMesh& mesh = meshPtr_();
    
    phiPtr_.reset(
        new surfaceScalarField(
            IOobject(
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("phi", dimensionSet(0, 3, -1, 0, 0, 0, 0), 0.0)));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void DASonicFoam::initSolver()
{
    /*
    Description:
        Initialize variables for DASolver
    */

    Info << "Initializing fields for DASonicFoam" << endl;
    Time& runTime = runTimePtr_();
    fvMesh& mesh = meshPtr_();
    argList& args = argsPtr_();
#include "createPimpleControlPython.H"
#include "createFieldsSonic.H"

    // read the RAS model from constant/turbulenceProperties
    const word turbModelName(
        IOdictionary(
            IOobject(
                "turbulenceProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false))
            .subDict("RAS")
            .lookup("RASModel"));
    daTurbulenceModelPtr_.reset(DATurbulenceModel::New(turbModelName, mesh, daOptionPtr_()));

#include "createAdjoint.H"


    // initialize fvSource and compute the source term
    const dictionary& allOptions = daOptionPtr_->getAllOptions();
    if (allOptions.subDict("fvSource").toc().size() != 0)
    {
        hasFvSource_ = 1;
        Info << "Initializing DASource" << endl;
        word sourceName = allOptions.subDict("fvSource").toc()[0];
        word fvSourceType = allOptions.subDict("fvSource").subDict(sourceName).getWord("type");
        daFvSourcePtr_.reset(DAFvSource::New(
            fvSourceType, mesh, daOptionPtr_(), daModelPtr_(), daIndexPtr_()));
    }

    // reduceIO does not write mesh, but if there is a FFD variable, set writeMesh to 1
    dictionary dvSubDict = daOptionPtr_->getAllOptions().subDict("inputInfo");
    forAll(dvSubDict.toc(), idxI)
    {
        word dvName = dvSubDict.toc()[idxI];
        if (dvSubDict.subDict(dvName).getWord("type") == "volCoord")
        {
            reduceIOWriteMesh_ = 1;
            break;
        }
    }
}

label DASonicFoam::solvePrimal()
{
    /*
    Description:
        Call the primal solver to get converged state variables
    */

    #include "createRefsSonic.H"

    daTurbulenceModelPtr_->updateIntermediateVariables();

    Info << "\nStarting time loop\n" << endl;

    label pimplePrintToScreen = 0;

    // reduce disk writes during adjoint run
    label reduceIO = daOptionPtr_->getAllOptions().subDict("unsteadyAdjoint").getLabel("reduceIO");
    wordList additionalOutput;
    if (reduceIO)
    {
        daOptionPtr_->getAllOptions().subDict("unsteadyAdjoint").readEntry<wordList>("additionalOutput", additionalOutput);
    }

    scalar endTime = runTime.endTime().value();
    scalar deltaT = runTime.deltaT().value();
    label nInstances = round(endTime / deltaT);

    label regModelFail = 0;
    label fail = 0;

    for (label iter = 1; iter <= nInstances; iter++)
    {
        ++runTime;

        printToScreen_ = this->isPrintTime(runTime, printIntervalUnsteady_);

        if (printToScreen_)
        {
            Info << "Time = " << runTime.timeName() << nl << endl;
            #include "compressibleCourantNo.H"
        }

        while (pimple.loop())
        {
            if (pimple.finalIter() && printToScreen_)
            {
                pimplePrintToScreen = 1;
            }
            else
            {
                pimplePrintToScreen = 0;
            }

            #include "UEqnSonic.H"
            #include "EEqnSonic.H"

            while (pimple.correct())
            {
                #include "pEqnSonic.H"
            }

            daTurbulenceModelPtr_->correct(pimplePrintToScreen);
            fail = daRegressionPtr_->compute();
        }

        regModelFail += fail;

        if (this->validateStates())
        {
            runTime.writeNow();
            mesh.write();
            return 1;
        }

        this->calcAllFunctions(printToScreen_);
        daRegressionPtr_->printInputInfo(printToScreen_);
        daTurbulenceModelPtr_->printYPlus(printToScreen_);
        this->printElapsedTime(runTime, printToScreen_);

        if (reduceIO && iter < nInstances)
        {
            this->writeAdjStates(reduceIOWriteMesh_, additionalOutput);
            daRegressionPtr_->writeFeatures();
        }
        else
        {
            runTime.write();
            daRegressionPtr_->writeFeatures();
        }
    }

    if (regModelFail != 0)
    {
        return 1;
    }

    primalFinalTimeIndex_ = runTime.timeIndex();
    mesh.write();

    Info << "End\n" << endl;

    return 0;
}
} // End namespace Foam

// ************************************************************************* //
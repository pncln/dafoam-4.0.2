#!/usr/bin/env python
"""
Run Python tests for DASonicPimpleFoam supersonic solver optimization integration
"""

from mpi4py import MPI
import os
import numpy as np
from testFuncs import *

import openmdao.api as om
from openmdao.api import Group
from mphys.multipoint import Multipoint
from dafoam.mphys.mphys_dafoam import DAFoamBuilderUnsteady
from mphys.scenario_aerodynamic import ScenarioAerodynamic
from pygeo.mphys import OM_DVGEOCOMP
from pygeo import geo_utils

gcomm = MPI.COMM_WORLD

os.chdir("./reg_test_files-main/ConvergentChannel")
if gcomm.rank == 0:
    os.system("rm -rf 0/* processor* *.bin")
    os.system("cp -r 0.compressible/* 0/")
    os.system("cp -r system.compressible.unsteady/* system/")
    os.system("cp -r constant/turbulenceProperties.sa constant/turbulenceProperties")
    os.system("cp -r constant/thermophysicalProperties.air constant/thermophysicalProperties")
    replace_text_in_file("system/fvSchemes", "meshWave;", "meshWaveFrozen;")

# supersonic aero setup
M0 = 2.0  # Mach number
T0 = 300.0  # Temperature [K]
p0 = 101325.0  # Pressure [Pa]
gamma = 1.4  # Heat capacity ratio
R = 287.0  # Gas constant [J/kg/K]

# Calculate supersonic flow conditions
c0 = np.sqrt(gamma * R * T0)  # Speed of sound
U0 = M0 * c0  # Velocity magnitude
rho0 = p0 / (R * T0)  # Density

daOptions = {
    "designSurfaces": ["walls"],
    "solverName": "DASonicPimpleFoam",
    "useAD": {"mode": "reverse", "seedIndex": 0, "dvName": "shape"},
    "primalBC": {
        "U0": {"variable": "U", "patches": ["inlet"], "value": [U0, 0.0, 0.0]},
        "p0": {"variable": "p", "patches": ["inlet"], "value": [p0]},
        "T0": {"variable": "T", "patches": ["inlet"], "value": [T0]},
        "rho0": {"variable": "rho", "patches": ["inlet"], "value": [rho0]},
        "useWallFunction": False,
    },
    "unsteadyAdjoint": {
        "mode": "timeAccurate",
        "PCMatPrecomputeInterval": 5,
        "PCMatUpdateInterval": 1,
        "readZeroFields": True,
        "additionalOutput": ["U", "p", "phi", "rho", "T"],
    },
    "shockCapturing": {
        "scheme": "artificialViscosity",
        "coefficient": 0.1,
        "useLocalTimeStep": 0,
    },
    "function": {
        "CD": {
            "type": "force",
            "source": "patchToFace",
            "patches": ["walls"],
            "directionMode": "fixedDirection",
            "direction": [1.0, 0.0, 0.0],
            "scale": 1.0 / (0.5 * rho0 * U0 * U0 * 1.0),
            "addToAdjoint": True,
        },
        "CL": {
            "type": "force",
            "source": "patchToFace",
            "patches": ["walls"],
            "directionMode": "fixedDirection",
            "direction": [0.0, 1.0, 0.0],
            "scale": 1.0 / (0.5 * rho0 * U0 * U0 * 1.0),
            "addToAdjoint": True,
        },
    },
    "adjoint": {
        "objFunc": "CD",
        "flowBoundaryConditions": {
            "bc0": {"patches": ["inlet"], "variable": "U", "value": [U0, 0.0, 0.0]},
            "bc1": {"patches": ["inlet"], "variable": "p", "value": [p0]},
            "bc2": {"patches": ["inlet"], "variable": "T", "value": [T0]},
            "bc3": {"patches": ["outlet"], "variable": "p", "value": [p0 * 0.5]},  # Supersonic outlet
        },
        "normalizeStates": {"U": U0, "p": p0, "rho": rho0, "T": T0, "phi": rho0 * U0},
        "referenceStates": {"U": U0, "p": p0, "rho": rho0, "T": T0, "phi": rho0 * U0},
        "adjEqnOption": {"gmresRelTol": 1.0e-6, "pcFillLevel": 1, "jacMatReOrdering": "rcm"},
    },
    "normalizeStates": {"U": U0, "p": p0, "rho": rho0, "T": T0, "phi": rho0 * U0},
    "referenceStates": {"U": U0, "p": p0, "rho": rho0, "T": T0, "phi": rho0 * U0},
    "primalMinResTol": 1e-8,
    "primalBC": {
        # Set supersonic boundary conditions
        "useWallFunction": False,
    },
    "fvSource": {},
    "inputInfo": {
        "omega": {"type": "design", "design": "shape", "patches": ["walls"], "size": 20},
    },
    "inputInfo": {
        "alpha": {"type": "parameter", "value": 0.0},
    },
    "multiPoint": False,
}

# mesh warping parameters
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openfoam",
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]],
}

# DVGeo
FFDFile = "./FFD/wingFFD.xyz"
DVGeo = geo_utils.createFFD(FFDFile, refAxis="y")

# select points
pts = DVGeo.getLocalIndex(0)
indexList = pts[:, 1].flatten()
PS = geo_utils.pointSet(pts, indexList, DVGeo)
DVGeo.addPointSet(PS, "wing")

daOptions["designVar"] = {}
nTwists = 6
for i in range(nTwists):
    daOptions["designVar"]["twist_%d" % i] = {"designVarType": "FFD"}

daOptions["designVar"]["alpha"] = {"designVarType": "AOA", "patches": ["inlet"], "flowAxis": "x", "normalAxis": "y"}

# DAFoam
DASolver = PYDAFOAM(options=daOptions, comm=gcomm)
DASolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
DASolver.printFamilyList()
DASolver.setMesh(mesh)
# set evalFuncs
evalFuncs = []
DASolver.setEvalFuncs(evalFuncs)

# DVConstraints
DVCon = geo_utils.DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface([DASolver.getTriangulatedMeshSurface(groupName=DASolver.designSurfacesGroup)])

# optFuncs
optFuncs.DASolver = DASolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm

if gcomm.rank == 0:
    os.system("rm -f 0.tar.gz && tar -czf 0.tar.gz 0 system")

time.sleep(1.0)

# Run
if __name__ == "__main__":
    DASolver.runColoring()
    xDV = DASolver.getVec()
    funcs = {}
    funcs, fail = optFuncs.calcObjCon(xDV)
    tol = 1.0e-4
    if gcomm.rank == 0:
        reg_write_dict(funcs, 1e-8, 1e-10)
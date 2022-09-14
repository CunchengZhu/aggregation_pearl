from operator import ge
import pymem3dg as dg
import pymem3dg.util as dg_util
import pymem3dg.broilerplate as dg_broil
import numpy as np
import copy
import pymem3dg.read as dg_read
from numpy.random import RandomState
from functools import partial


def initialConditionsByMatrices():
    lengthScale = 1
    face, vertex = dg.getCylinder(radius = lengthScale, radialSubdivision = 16, axialSubdivision = 60)
    proteinDensity = np.ones(np.shape(vertex)[0])
    velocity = np.zeros(np.shape(vertex))

    initialConditions = {
        "topologyMatrix": face,
        "vertexMatrix": vertex,
        "referenceVertexMatrix": vertex,
        "proteinDensity": proteinDensity,
        "velocity": velocity,
    }
    return initialConditions, lengthScale


def continuationByNc():
    """trajFile construction"""
    outputDir = "."
    trajFile = outputDir + "//traj.nc"
    FRAME = dg_read.sizeOf(trajFile) - 1
    initialConditions = {"trajFile": trajFile, "startingFrame": FRAME}
    return initialConditions


def parameters():
    p = dg.Parameters()
    p.point.index = 0

    p.boundary.shapeBoundaryCondition = "roller"
    p.boundary.proteinBoundaryCondition = "none"

    p.variation.isProteinVariation = False
    p.variation.isProteinConservation = False
    p.variation.isShapeVariation = True

    p.bending.Kd = 0
    p.bending.Kdc = 0
    p.bending.Kb = 0
    p.bending.Kbc = 0.1
    p.bending.H0c = 1

    p.tension.A_res = 0
    p.tension.setForm(partial(dg_broil.constantSurfaceTensionModel, tension = 1e-7))

    p.osmotic.V_res = 4/ 3 * 3.14
    p.osmotic.setForm(partial(dg_broil.preferredVolumeOsmoticPressureModel, strength = 100, preferredVolume = 15.6736 * 0.5))
    return p


if __name__ == "__main__":
    p = parameters()

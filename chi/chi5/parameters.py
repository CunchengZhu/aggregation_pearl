import pymem3dg as dg
import pymem3dg.util as dg_util
import numpy as np
import copy
import pymem3dg.read as dg_read
from numpy.random import RandomState

def scalingVariables():
    # temp variable
    xi = 1
    A_bar = 12.4866
    R_bar = np.sqrt(A_bar / 4 / np.pi)
    Kb = 8.22e-5
    h = 1e-7 * (xi * R_bar**2 / Kb)
    return xi, A_bar, R_bar, Kb, h


def periodic(vertexPositions, vertexDualAreas, time, geodesicDistance):
    freq = 12
    totalHeight = np.max(vertexPositions[:, 2]) - np.min(vertexPositions[:, 2])

    direction = copy.deepcopy(vertexPositions)
    direction[:, 2] = 0
    direction = dg_util.rowwiseNormalize(direction)

    xi, A_bar, R_bar, Kb, h = scalingVariables()
    Kf = 0.05 * (Kb / R_bar)

    magnitude = Kf * (
        1 + np.sin(freq * 2 * np.pi / totalHeight * vertexPositions[:, 2])
    )

    return dg_util.rowwiseScaling(magnitude, direction)


def point(vertexPositions, vertexDualAreas, time, geodesicDistances):
    xi, A_bar, R_bar, Kb, h = scalingVariables()
    decayTime = 100000 * h
    std = 0.02 * R_bar
    Kf = 5 * (Kb / R_bar)
    direction = dg_util.rowwiseNormalize(vertexPositions)
    magnitude = (
        Kf
        * np.exp(-time / decayTime)
        * dg_util.gaussianDistribution(geodesicDistances, 0, std)
    )

    return dg_util.rowwiseScaling(magnitude, direction)


def initialConditionsByMatrices():
    """matrix construction"""
    # face, vertex = dg.getIcosphere(1, 3)
    face, vertex = dg_read.readMeshByPly("../inputMesh.ply")
    # vertex = dg_util.sphericalHarmonicsPerturbation(vertex, 5, 6, 0.1)
    # refVertex = dg_util.sphericalHarmonicsPerturbation(vertex, 5, 6, 0.1)
    refVertex = vertex
    prng = RandomState(1234567890)
    proteinDensity = prng.rand(np.shape(vertex)[0])
    velocity = np.zeros(np.shape(vertex))
    FRAME = 0
    return face, vertex, refVertex, proteinDensity, velocity, FRAME


def continuationByNc():
    """trajFile construction"""
    outputDir = "."
    trajFile = outputDir + "//traj.nc"
    FRAME = dg_read.sizeOf(trajFile) - 1
    return trajFile, FRAME


def parameters(xi, A_bar, R_bar, Kb):
    p = dg.Parameters()

    p.proteinMobility = 0.01 * (1 / xi / R_bar**2)
    p.temperature = 0

    p.point.pt = [0, 0, 10 * R_bar]
    p.point.isFloatVertex = False

    p.protein.profile = "none"
    p.protein.geodesicProteinDensityDistribution = [-1]
    p.protein.proteinInteriorPenalty = 0

    p.boundary.shapeBoundaryCondition = "none"
    p.boundary.proteinBoundaryCondition = "none"

    p.variation.isProteinVariation = True
    p.variation.isProteinConservation = True
    p.variation.isShapeVariation = True
    p.variation.geodesicMask = -1

    p.bending.Kd = 0
    p.bending.Kdc = 0
    p.bending.Kb = Kb
    p.bending.Kbc = 0
    p.bending.H0c = 12 / R_bar

    p.tension.isConstantSurfaceTension = False
    p.tension.Ksg = 12000 * (Kb / R_bar**2)
    p.tension.A_res = 0
    p.tension.At = A_bar
    p.tension.lambdaSG = 0

    p.adsorption.epsilon = 0 * Kb / R_bar**2

    p.aggregation.chi = 5 * Kb / R_bar**2

    p.entropy.xi = 1 * Kb / R_bar**2

    p.osmotic.isPreferredVolume = True
    p.osmotic.isConstantOsmoticPressure = False
    p.osmotic.Kv = 500 * Kb
    p.osmotic.V_res = 0
    p.osmotic.n = 1
    p.osmotic.Vt = 0.85 * (4 / 3 * np.pi * R_bar**3)
    p.osmotic.cam = -1
    p.osmotic.lambdaV = 0

    p.dirichlet.eta = 0.1 * Kb

    p.selfAvoidance.d = 0.001
    p.selfAvoidance.mu = 0
    p.selfAvoidance.p = 0
    p.selfAvoidance.n = 2

    p.dpd.gamma = 0

    p.external.setForm(None)

    p.spring.Kst = 0
    p.spring.Ksl = 0
    p.spring.Kse = 0
    return p


if __name__ == "__main__":
    p = parameters()

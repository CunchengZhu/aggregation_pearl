from operator import ge
import pymem3dg as dg
import pymem3dg.util as dg_util
import numpy as np
import copy
import pymem3dg.read as dg_read
from numpy.random import RandomState


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
    _, lengthScale = initialConditionsByMatrices()
    decayTime = 30
    std = 0.05 * lengthScale
    Kf = 0.2
    direction = dg_util.rowwiseNormalize(vertexPositions)
    magnitude = (
        Kf
        * np.exp(-time / decayTime)
        * dg_util.gaussianDistribution(geodesicDistances, 0, std)
    )

    return dg_util.rowwiseScaling(magnitude, direction)


def initialConditionsByMatrices():
    lengthScale = 6
    face, vertex = dg.getIcosphere(lengthScale, 3)
    # face, vertex = dg_read.readMeshByPly("../inputMesh.ply")
    # vertex = dg_util.sphericalHarmonicsPerturbation(vertex, 5, 6, 0.1 * lengthScale)
    # refVertex = dg_util.sphericalHarmonicsPerturbation(vertex, 5, 6, 0.1)

    # prng = RandomState(1234567890)
    # proteinDensity = prng.rand(np.shape(vertex)[0])
    # proteinDensity = np.ones(np.shape(vertex)[0]) * 0.2
    p = parameters()
    geodesicDistance = dg_util.getGeodesicDistance(face, vertex, p.point.pt)
    proteinDensity = dg_util.tanhDistribution(
        geodesicDistance,
        p.protein.tanhSharpness,
        p.protein.geodesicProteinDensityDistribution[0],
    )
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

    p.proteinMobility = 0
    p.temperature = 0

    p.point.pt = [0]
    p.point.isFloatVertex = False

    p.protein.profile = "tanh"
    p.protein.geodesicProteinDensityDistribution = [2, 2, 1, 0]
    p.protein.proteinInteriorPenalty = 0
    p.protein.tanhSharpness = 10

    p.boundary.shapeBoundaryCondition = "none"
    p.boundary.proteinBoundaryCondition = "none"

    p.variation.isProteinVariation = False
    p.variation.isProteinConservation = False
    p.variation.isShapeVariation = True
    p.variation.geodesicMask = -1

    p.bending.Kd = 0
    p.bending.Kdc = 0.1
    p.bending.Kb = 0.1
    p.bending.Kbc = 0.2
    p.bending.H0c = 3

    p.tension.isConstantSurfaceTension = False
    p.tension.Ksg = 1e4
    p.tension.A_res = 0
    p.tension.At = 453.874  # 4 * np.pi * R**2
    p.tension.lambdaSG = 0

    p.adsorption.epsilon = 0

    p.aggregation.chi = 2 * 0

    p.entropy.xi = 1 * 0

    p.osmotic.isPreferredVolume = True
    p.osmotic.isConstantOsmoticPressure = False
    p.osmotic.Kv = 500
    p.osmotic.V_res = 0
    p.osmotic.n = 1
    p.osmotic.Vt = 0.9 * 897  # (4 / 3 * np.pi * R**3)
    p.osmotic.cam = -1
    p.osmotic.lambdaV = 0

    p.dirichlet.eta = 2 * 0

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

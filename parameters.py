import pymem3dg as dg
import pymem3dg.util as dg_util
import numpy as np
import copy


def scalingVariables():
    # temp variable
    xi = 1
    kT = 0.411e-5
    R_bar = 0.05  # connector radius
    Kb = 50 * kT
    Ksg = 0.01
    h = xi * R_bar ** 2 / Kb
    return xi, Ksg, Kb, R_bar, h


def periodic(vertexPositions, vertexDualAreas, time, geodesicDistance):
    freq = 12
    totalHeight = np.max(vertexPositions[:, 2]) - np.min(vertexPositions[:, 2])

    direction = copy.deepcopy(vertexPositions)
    direction[:, 2] = 0
    direction = dg_util.rowwiseNormalize(direction)

    xi, Ksg, Kb, R_bar, h = scalingVariables()
    Kf = 0.05 * (Kb / R_bar)

    magnitude = Kf * (
        1 + np.sin(freq * 2 * np.pi / totalHeight * vertexPositions[:, 2])
    )

    return dg_util.rowwiseScaling(magnitude, direction)


def point(vertexPositions, vertexDualAreas, time, geodesicDistances):
    xi, Ksg, Kb, R_bar, h = scalingVariables()
    decayTime = 100 * h
    std = 1 * R_bar
    Kf = 1 * (Kb / R_bar)

    direction = -dg_util.rowwiseNormalize(vertexPositions - [0, 0, 10 * R_bar])
    magnitude = (
        Kf
        * np.exp(-time / decayTime)
        * dg_util.gaussianDistribution(geodesicDistances, 0, std)
    )

    return dg_util.rowwiseScaling(magnitude, direction)


def parameters(xi, R_bar, Kb, Ksg):
    p = dg.Parameters()

    p.proteinMobility = 0 * (1 / xi / R_bar ** 2)
    p.temperature = 0

    p.point.pt = [R_bar, R_bar, 10 * R_bar]
    p.point.isFloatVertex = False

    p.protein.profile = "none"
    p.protein.geodesicProteinDensityDistribution = [-1]
    p.protein.proteinInteriorPenalty = 0

    p.boundary.shapeBoundaryCondition = "fixed"
    p.boundary.proteinBoundaryCondition = "none"

    p.variation.isProteinVariation = False
    p.variation.isShapeVariation = True
    p.variation.geodesicMask = -1

    p.bending.Kd = 0
    p.bending.Kdc = 0
    p.bending.Kb = 0
    p.bending.Kbc = 1 * Kb
    p.bending.H0c = ((2 * R_bar) ** (-2) - Ksg / Kb) ** 0.5

    p.tension.isConstantSurfaceTension = True
    p.tension.Ksg = Ksg
    p.tension.A_res = 0
    p.tension.At = -1
    p.tension.lambdaSG = 0

    p.adsorption.epsilon = 0

    p.aggregation.chi = 0

    p.osmotic.isPreferredVolume = True
    p.osmotic.isConstantOsmoticPressure = False
    p.osmotic.Kv = 500 * Kb
    p.osmotic.V_res = 0 * (4 / 3 * np.pi * R_bar ** 3)
    p.osmotic.n = 1
    p.osmotic.Vt = 90 * (4 / 3 * np.pi * R_bar ** 3)
    p.osmotic.cam = -1
    p.osmotic.lambdaV = 0

    p.dirichlet.eta = 0

    p.selfAvoidance.d = 0.001
    p.selfAvoidance.mu = 0
    p.selfAvoidance.p = 0
    p.selfAvoidance.n = 2

    p.dpd.gamma = 0

    p.external.setForm(point)
    return p


if __name__ == "__main__":
    p = parameters()

from inspect import Parameter
from operator import ge
import pymem3dg as dg
import pymem3dg.util as dg_util
import numpy as np
import copy
import pymem3dg.read as dg_read
from numpy.random import RandomState
from functools import partial
import numpy.typing as npt


lengthScale = 6

def initialConditionsByMatrices():
    face, vertex = dg.getIcosphere(lengthScale, 3)
    notableVertex = 0
    geometry = dg.Geometry(face, vertex, vertex, notableVertex)
    geodesicDistance = geometry.computeGeodesicDistance()
    p = parameters()
    proteinDensity = p.protein.form(
        0, np.zeros([np.shape(vertex)[0], 1]), geodesicDistance
    )  # dummy placeholder args time & mean curvature
    velocity = np.zeros(np.shape(vertex))

    initialConditions = {
        "geometry": geometry,
        "proteinDensity": proteinDensity,
        "velocity": velocity,
        "parameters": p,
    }
    return initialConditions

def meshProcessorSetting(meshProcessor: dg.MeshProcessor):
    meshProcessor.meshMutator.isShiftVertex = True
    meshProcessor.meshMutator.flipNonDelaunay = True
    # meshProcessor.meshMutator.splitLarge = True
    meshProcessor.meshMutator.splitFat = True
    meshProcessor.meshMutator.splitSkinnyDelaunay = True
    meshProcessor.meshMutator.splitCurved = True
    meshProcessor.meshMutator.minimumEdgeLength = 0.02 * lengthScale
    meshProcessor.meshMutator.maximumEdgeLength = 0.2 * lengthScale
    meshProcessor.meshMutator.curvTol = 0.1 / lengthScale
    meshProcessor.meshMutator.collapseSkinny = True
    meshProcessor.meshMutator.collapseSmall = True
    meshProcessor.meshMutator.collapseFlat = True
    # meshProcessor.meshMutator.targetFaceArea = 0.0003 * lengthScale **2
    meshProcessor.meshMutator.isSmoothenMesh = True
    return meshProcessor

def continuationByNc(frame: int = None):
    trajFile = ".//traj.nc"
    if frame == None:
        frame = dg_read.sizeOf(trajFile) - 1

    p = parameters()
    geometry = dg.Geometry(trajFile, frame)
    initialConditions = {
        "geometry": geometry,
        "trajFile": trajFile,
        "startingFrame": frame,
        "parameters": p,
    }
    return initialConditions


def prescribeGeodesicPoteinDensityDistribution(
    time: float,
    vertexMeanCuravtures: npt.NDArray[np.float64],
    geodesicDistance: npt.NDArray[np.float64],
    sharpness: float,
    radius: float,
):
    """form function that prescribe geodesic protein density profile on mesh

    Args:
        time (float): time of the system
        vertexMeanCuravtures (npt.NDArray[np.float64]): vertex mean curvature
        geodesicDistance (npt.NDArray[np.float64]): vertex geodesic curvature
        sharpness (float): sharpness of transition used in tanh function
        radius (float): center of transition of the tanh function

    Returns:
        _type_: _description_
    """
    return dg_util.tanhDistribution(
        x=geodesicDistance, sharpness=sharpness, center=radius
    )


def preferredAreaSurfaceTensionModel(
    area: float, modulus: float, preferredArea: float, reservoirArea: float = 0
):
    """harmonic potential type model with a preferred area

    Args:
        area (float): total surface area
        modulus (float): stretching modulus
        preferredArea (float): value of preferred area

    Returns:
        tuple: surface tension and surface energy
    """
    area_difference = (area + reservoirArea) - preferredArea
    tension = modulus * area_difference / preferredArea
    energy = tension * area_difference / 2
    return (tension, energy)


def preferredVolumeOsmoticPressureModel(
    volume: float, strength: float, preferredVolume: float, reservoirVolume: float = 0
):
    """harmonic potential type model with a preferred volume

    Args:
        volume (float): enclosed volume
        strength (float): osmotic stregth coefficient 
        preferredVolume (float): value of preferred volume

    Returns:
        tuple: osmotic pressure and pressure energy
    """
    volume_difference = (volume + reservoirVolume) - preferredVolume
    pressure = -(strength * volume_difference / preferredVolume / preferredVolume)
    return (pressure, -pressure * volume_difference / 2)


def parameters():
    p = dg.Parameters()

    p.proteinMobility = 0
    p.temperature = 0

    p.protein.proteinInteriorPenalty = 0
    p.protein.form = partial(
        prescribeGeodesicPoteinDensityDistribution, sharpness=10, radius=1
    )

    p.boundary.shapeBoundaryCondition = "none"
    p.boundary.proteinBoundaryCondition = "none"

    p.variation.isProteinVariation = False
    p.variation.isProteinConservation = False
    p.variation.isShapeVariation = True
    p.variation.geodesicMask = -1

    p.bending.Kd = 0.005
    p.bending.Kdc = 0.02
    p.bending.Kb = 0.1
    p.bending.Kbc = 0.2
    p.bending.H0c = 2

    p.tension.form = partial(
        preferredAreaSurfaceTensionModel,
        modulus=1e4,
        preferredArea=453.874,
        reservoirArea=0,
    )

    p.adsorption.epsilon = 0

    p.aggregation.chi = 2 * 0

    p.entropy.xi = 1 * 0

    p.osmotic.form = partial(
        preferredVolumeOsmoticPressureModel,
        strength=3000,
        preferredVolume=0.9 * 909.236,  # (4 / 3 * np.pi * R**3)
        reservoirVolume=0,
    )

    p.dirichlet.eta = 2 * 0

    p.selfAvoidance.d = 0.001
    p.selfAvoidance.mu = 0
    p.selfAvoidance.p = 0
    p.selfAvoidance.n = 2

    p.dpd.gamma = 0

    p.external.form = None

    p.spring.Kst = 0
    p.spring.Ksl = 0
    p.spring.Kse = 0
    return p


if __name__ == "__main__":
    p = parameters()

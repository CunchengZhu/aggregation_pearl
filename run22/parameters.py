from inspect import Parameter
from operator import ge
import pymem3dg as dg
import pymem3dg.util as dg_util
import numpy as np
import copy
import pymem3dg.read.netcdf as dg_nc
from numpy.random import RandomState
from functools import partial
import numpy.typing as npt


lengthScale = 1

def initialConditionsByMatrices():
    face, vertex = dg.getIcosphere(lengthScale, 3)
    notableVertex = [False] * np.shape(vertex)[0]
    notableVertex[0] = True
    notableVertex[100] = True
    geometry = dg.Geometry(face, vertex, vertex, notableVertex)
    proteinDensity = dg_util.tanhDistribution(
        x=geometry.computeGeodesicDistance(), sharpness=20, center=0.3
    )
    p = parameters()
    # proteinDensity = np.ones(np.shape(vertex)[0]) * 0.1
    velocity = np.zeros(np.shape(vertex))

    initialConditions = {
        "geometry": geometry,
        "proteinDensity": proteinDensity,
        "velocity": velocity,
        "parameters": p,
    }
    return initialConditions

def meshProcessorSetting(meshProcessor: dg.MeshProcessor):
    meshProcessor.meshMutator.mutateMeshPeriod = 10
    meshProcessor.meshMutator.isShiftVertex = True
    meshProcessor.meshMutator.flipNonDelaunay = True
    # meshProcessor.meshMutator.splitLarge = True
    meshProcessor.meshMutator.splitFat = True
    meshProcessor.meshMutator.splitSkinnyDelaunay = True
    meshProcessor.meshMutator.splitCurved = True
    meshProcessor.meshMutator.minimumEdgeLength = 0.02 * lengthScale
    meshProcessor.meshMutator.maximumEdgeLength = 0.2 * lengthScale
    meshProcessor.meshMutator.curvTol = 0.004 / lengthScale
    meshProcessor.meshMutator.collapseSkinny = True
    meshProcessor.meshMutator.collapseSmall = True
    meshProcessor.meshMutator.collapseFlat = True
    # meshProcessor.meshMutator.targetFaceArea = 0.0003 * lengthScale **2
    meshProcessor.meshMutator.isSmoothenMesh = True
    return meshProcessor

def continuationByNc(frame: int = None):
    trajFile = ".//traj.nc"
    if frame == None:
        frame = dg_nc.sizeOf(trajFile) - 1

    p = parameters()
    geometry = dg.Geometry(trajFile, frame)
    initialConditions = {
        "geometry": geometry,
        "trajFile": trajFile,
        "startingFrame": frame,
        "parameters": p,
    }
    return initialConditions


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

    p.proteinMobility = 0.01
    p.temperature = 0

    p.variation.isProteinVariation = True
    p.variation.isProteinConservation = True
    p.variation.isShapeVariation = True
    
    p.bending.Kd = 2e-7
    p.bending.Kb = 8.22e-5
    p.bending.Kbc = 8.22e-5
    p.bending.H0c = 10

    p.tension.form = partial(
        preferredAreaSurfaceTensionModel,
        modulus=1,
        preferredArea=12.5025,
        reservoirArea=0,
    )

    p.adsorption.epsilon = -1e-3

    p.osmotic.form = partial(
        preferredVolumeOsmoticPressureModel,
        strength=0.05,
        preferredVolume=0.7 * 4.15889,  # (4 / 3 * np.pi * R**3)
        reservoirVolume=0,
    )

    p.dirichlet.eta = 3 * p.bending.Kb
    
    p.aggregation.chi = 10 * p.bending.Kb

    return p


if __name__ == "__main__":
    p = parameters()

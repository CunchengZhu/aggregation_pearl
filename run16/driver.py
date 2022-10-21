import pymem3dg as dg
import pymem3dg.util as dg_util
import pymem3dg.visual as dg_vis
import pymem3dg.read as dg_read
import numpy as np
import parameters as ps

def localRun():
    # CONTINUE = True
    CONTINUE = False
    ####################################################
    #                 Initialize pathes                #
    ####################################################
    outputDir = "."
    ####################################################
    #            Arguments for System                    #
    ####################################################
    parameters = ps.parameters()
    if CONTINUE:
        initialConditions, lengthScale = ps.initialConditionsByMatrices()
        initialConditions = ps.continuationByNc()
    else:
        initialConditions, lengthScale = ps.initialConditionsByMatrices()
    ####################################################
    #                 System                           #
    ####################################################
    """ System construction """
    arguments = initialConditions
    arguments["parameters"] = parameters
    g = dg.System(**arguments)

    """ Mesh processor """
    g.meshProcessor.meshMutator.isShiftVertex = True
    g.meshProcessor.meshMutator.flipNonDelaunay = True
    # g.meshProcessor.meshMutator.splitLarge = True
    g.meshProcessor.meshMutator.splitFat = True
    g.meshProcessor.meshMutator.splitSkinnyDelaunay = True
    g.meshProcessor.meshMutator.splitCurved = True
    g.meshProcessor.meshMutator.minimumEdgeLength = 0.001 * lengthScale
    g.meshProcessor.meshMutator.curvTol = 0.1 / lengthScale
    g.meshProcessor.meshMutator.collapseSkinny = True
    g.meshProcessor.meshMutator.collapseSmall = True
    # g.meshProcessor.meshMutator.collapseFlat = True
    # g.meshProcessor.meshMutator.targetFaceArea = 0.0003 * lengthScale **2
    g.meshProcessor.meshMutator.isSmoothenMesh = True
    """ System initialization """
    g.initialize(nMutation=0, ifMute=False)
    # g.testForceComputation(h)
    ####################################################
    #          Time integration / Optimization         #
    ####################################################
    """ Integrator construction """
    if CONTINUE:
        frame = initialConditions["startingFrame"]
    else:
        frame = 0
    fe = dg.Euler(
        system=g,
        characteristicTimeStep=5e-4,
        totalTime=3000000,
        savePeriod=500 * 5e-4,
        tolerance=0,
        outputDirectory=outputDir,
        frame=frame,
    )
    """ settings """
    fe.updateGeodesicsPeriod = 10
    fe.processMeshPeriod = 10
    # fe.fluctuatePeriod = 10
    # fe.fluctuateAmplitude = 0.001
    fe.isBacktrack = True
    # fe.ifAdaptiveStep = False
    """ Verbosity """
    fe.ifPrintToConsole = True
    fe.ifOutputTrajFile = True
    # fe.ifOutputMeshFile = True

    fe.integrate()


if __name__ == "__main__":
    localRun()

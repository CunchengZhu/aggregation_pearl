import pymem3dg as dg
import pymem3dg.util as dg_util
import pymem3dg.visual as dg_vis
import pymem3dg.read as dg_read
import numpy as np
import parameters as ps

# CONTINUE = True
CONTINUE = False


def localRun(parameters, CONTINUE=True):
    ####################################################
    #                 Initialize pathes                #
    ####################################################
    outputDir = "."
    ####################################################
    #                 Parameters                       #
    ####################################################
    xi, A_bar, R_bar, Kb, h = ps.scalingVariables()
    ####################################################
    #            Initial conditions                    #
    ####################################################
    if CONTINUE:
        trajFile, FRAME = ps.continuationByNc()
    else:
        (
            face,
            vertex,
            refVertex,
            proteinDensity,
            velocity,
            FRAME,
        ) = ps.initialConditionsByMatrices()
    ####################################################
    #                 System                           #
    ####################################################
    """ System construction """
    if CONTINUE:
        g = dg.System(trajFile, FRAME, parameters)
    else:
        g = dg.System(face, vertex, refVertex, proteinDensity, velocity, parameters)
    """ Mesh processor """
    g.meshProcessor.meshMutator.isShiftVertex = True
    g.meshProcessor.meshMutator.flipNonDelaunay = True
    # g.meshProcessor.meshMutator.splitLarge = True
    g.meshProcessor.meshMutator.splitFat = True
    g.meshProcessor.meshMutator.splitSkinnyDelaunay = True
    g.meshProcessor.meshMutator.splitCurved = True
    g.meshProcessor.meshMutator.minimumEdgeLength = 0.001 * R_bar
    g.meshProcessor.meshMutator.curvTol = 0.004 / R_bar
    g.meshProcessor.meshMutator.collapseSkinny = True
    g.meshProcessor.meshMutator.collapseSmall = True
    g.meshProcessor.meshMutator.collapseFlat = True
    g.meshProcessor.meshMutator.targetFaceArea = 0.0003 * R_bar**2
    g.meshProcessor.meshMutator.isSmoothenMesh = True
    """ System initialization """
    g.initialize(nMutation=0, ifMute=False)
    # g.testForceComputation(h)
    ####################################################
    #          Time integration / Optimization         #
    ####################################################
    """ Integrator construction """
    fe = dg.Euler(
        system=g,
        characteristicTimeStep= 3 * h,
        totalTime=3000000 * h,
        savePeriod=1000 * h,
        tolerance=0.1 * (Kb / R_bar),
        outputDirectory=outputDir,
        frame=FRAME,
    )
    """ settings """
    fe.updateGeodesicsPeriod = 100
    fe.processMeshPeriod = 10000
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
    xi, A_bar, R_bar, Kb, h = ps.scalingVariables()
    parameters = ps.parameters(xi, A_bar, R_bar, Kb)
    localRun(parameters, CONTINUE=CONTINUE)

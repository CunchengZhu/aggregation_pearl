import pymem3dg as dg
import pymem3dg.util as dg_util
import pymem3dg.visual as dg_vis
import pymem3dg.read as dg_read
import numpy as np
import parameters as ps

def localRun():
    CONTINUE = True
    ####################################################
    #                 Initialize pathes                #
    ####################################################
    outputDir = "."
    ####################################################
    #            Arguments for System                    #
    ####################################################
    if CONTINUE:
        initialConditions = ps.continuationByNc()
    else:
        initialConditions = ps.initialConditionsByMatrices()
    ####################################################
    #                 System                           #
    ####################################################
    """ System construction """
    g = dg.System(**initialConditions)
    g.meshProcessor = ps.meshProcessorSetting(g.meshProcessor)
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

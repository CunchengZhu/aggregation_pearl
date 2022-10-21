import pymem3dg as dg
import parameters as ps


def localRun():
    CONTINUE = False
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
    g.initialize(ifMutateMesh=False, ifMute=False)

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
        characteristicTimeStep=0.01,
        totalTime=50000,
        savePeriod=1000 * 0.25,
        tolerance=0,
        outputDirectory=outputDir,
        frame=frame,
    )
    """ settings """
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

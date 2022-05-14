import pymem3dg as dg
import pymem3dg.util as dg_util
import pymem3dg.visual as dg_vis
import pymem3dg.read as dg_read
import numpy as np
import parameters

####################################################
#                 Initialize pathes                #
####################################################
outputDir = "."
####################################################
#                 Parameters                       #
####################################################
xi, A_bar, R_bar, Kb, h = parameters.scalingVariables()
p = parameters.parameters(xi, A_bar, R_bar, Kb)
####################################################
#            Initial conditions                    #
####################################################
(
    face,
    vertex,
    proteinDensity,
    velocity,
    FRAME
) = parameters.initialConditionsByMatrices()
# trajFile, FRAME = parameters.continuationByNc()
####################################################
#                 System                           #
####################################################
""" System construction """
g = dg.System(face, vertex, proteinDensity, velocity, p)
# g = dg.System(trajFile, FRAME, p)
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
    characteristicTimeStep=h,
    totalTime=1000000 * h,
    savePeriod=200 * h,
    tolerance=0.1 * (Kb / R_bar),
    outputDirectory=outputDir,
    frame=FRAME,
)
""" settings """
fe.updateGeodesicsPeriod = 100
fe.processMeshPeriod = 20
# fe.fluctuatePeriod = 10
# fe.fluctuateAmplitude = 0.001
fe.isBacktrack = True
# fe.ifAdaptiveStep = False
""" Verbosity """
fe.ifPrintToConsole = True
fe.ifOutputTrajFile = True
# fe.ifOutputMeshFile = True

fe.integrate()
####################################################
#                       Visualization              #
####################################################
# dg_vis.animate(outputDir+"/traj.nc", meanCurvature = True)
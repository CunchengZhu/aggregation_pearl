import pymem3dg as dg
import pymem3dg.util as dg_util
import pymem3dg.visual as dg_vis
import numpy as np
import parameters

####################################################
#                 Initialize pathes                #
####################################################
outputDir = "."
####################################################
#                 Parameters                       #
####################################################
""" Import from file """
xi, Ksg, Kb, R_bar, h = parameters.scalingVariables()
p = parameters.parameters(xi, R_bar, Kb, Ksg)
####################################################
#            Initial conditions                    #
####################################################
""" Built-in construction """
print("radius is ", 0.5 * (p.tension.Ksg / p.bending.Kbc + p.bending.H0c ** 2) ** (-0.5))
face, vertex = dg.getCylinder(
    radius=0.5 * (p.tension.Ksg / p.bending.Kbc + p.bending.H0c ** 2) ** (-0.5),
    radialSubdivision=int(
        (p.tension.Ksg / p.bending.Kbc + p.bending.H0c ** 2) ** (-0.5)
        / (p.tension.Ksg / Kb + p.bending.H0c ** 2) ** (-0.5)
        * 10
    ),
    axialSubdivision=60,
    frequency=7.5,
    amplitude=0,
)
""" input construction """
trajFile = outputDir + "//traj.nc"
# inputMesh = outputDir + "/temp7/frame1780.ply"
""" additional initial condition """
proteinDensity = np.ones(np.shape(vertex)[0]) * 0.5
velocity = np.zeros(np.shape(vertex))
####################################################
#                 System                           #
####################################################
FRAME = 0
""" System construction """
g = dg.System(face, vertex, p)
# g = dg.System(Face, Vertex, proteinDensity, velocity, p)
# g = dg.System(trajFile, FRAME, p)
""" Mesh processor """
g.meshProcessor.meshMutator.isShiftVertex = True
g.meshProcessor.meshMutator.flipNonDelaunay = True
# g.meshProcessor.meshMutator.splitLarge = True
# g.meshProcessor.meshMutator.splitFat = True
g.meshProcessor.meshMutator.splitSkinnyDelaunay = True
g.meshProcessor.meshMutator.splitCurved = True
g.meshProcessor.meshMutator.minimumEdgeLength = 0.001 * R_bar
g.meshProcessor.meshMutator.curvTol = 0.0005 / R_bar
g.meshProcessor.meshMutator.collapseSkinny = True
# g.meshProcessor.meshMutator.collapseSmall = True
g.meshProcessor.meshMutator.collapseFlat = True
g.meshProcessor.meshMutator.targetFaceArea = 0.01 * R_bar ** 2
g.meshProcessor.meshMutator.isSmoothenMesh = True
# g.meshProcessor.meshRegularizer.Kst = 0.1 # 2e-6
# g.meshProcessor.meshRegularizer.Ksl = 0
# g.meshProcessor.meshRegularizer.Kse = 0
# g.meshProcessor.meshRegularizer.readReferenceData(icoFace, icoVertex, 0)
""" System initialization """
g.initialize(nMutation=0, ifMute=False)
####################################################
#          Time integration / Optimization         #
####################################################
""" Integrator construction """
fe = dg.Euler(
    system=g,
    characteristicTimeStep=1e-2 * h,
    totalTime=1000 * h,
    savePeriod=10 * h,
    tolerance=0.1 * (Kb / R_bar),
    outputDirectory=outputDir,
    frame=FRAME,
)
""" settings """
fe.updateGeodesicsPeriod = 20
fe.processMeshPeriod = 20
# fe.fluctuatePeriod = 10
# fe.fluctuateAmplitude = 0.001
fe.isBacktrack = True
fe.ifAdaptiveStep = False
""" Verbosity """
fe.ifPrintToConsole = True
fe.ifOutputTrajFile = True
# fe.ifOutputMeshFile = True

fe.integrate()
####################################################
#                       Visualization              #
####################################################
# dg_vis.animate(outputDir+"/traj.nc", meanCurvature = True)

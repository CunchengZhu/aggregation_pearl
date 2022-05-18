import pymem3dg.visual as dg_vis
import imp

####################################################
#                 Initialize pathes                #
####################################################
# folder = "/beads_on_pearl"
# folder = "/intermediary"
# folder = r"H:/Shared drives/Rangamani Lab Drive/Cuncheng Zhu/manuscript_in_process/feng_pearling/repo/aggregation_pearl/chi20"
folder = "../run5"
var = 8
subFolder = folder + f"/chi{var}"
# folder = "/no_pressure"
# folder = "/tube"
parameterFile = imp.load_source("module.name", folder + "/parameters.py")
xi, Ksg, Kb, R_bar, h = parameterFile.scalingVariables()
parameters = parameterFile.parameters(xi, R_bar, Kb, Ksg)
parameters.aggregation.chi = var * Kb / R_bar**2
trajNc = subFolder + "/traj.nc"
####################################################
#                  visualize .ply                  #
####################################################
# dg_vis.visualizePly(plyFile, meanCurvature=True, gaussianCurvature=True)

####################################################
#            visualize animated .nc                #
####################################################
dg_vis.animate(
    trajNc,
    parameters=parameters,
    meanCurvature=True,
    gaussianCurvature=True,
    bendingForce=True,
    externalForce=True,
    mechanicalForce=True,
    capillaryForce=True,
    lineCapillaryForce=True,
    osmoticForce=True,
    springForce=True,
    entropyForce=True,
    chemicalPotential=True,
    bendingPotential=True,
    diffusionPotential=True,
    aggregationPotential=True,
    adsorptionPotential=True,
    entropyPotential=True
)

####################################################
#            visualize energy                      #
####################################################
dg_vis.plotEnergy(
    trajNc,
    parameters,
    potentialEnergy=True,
    # kineticEnergy=True,
    # totalEnergy=True,
    bendingEnergy=True,
    # externalWork=True,
    # deviatoricEnergy=True,
    surfaceEnergy=True,
    pressureEnergy=True,
    # adsorptionEnergy=True,
    aggregationEnergy=True,
    entropyEnergy=True,
    # edgeSpringEnergy=True,
    # faceSpringEnergy=True,
    # lcrSpringEnergy=True,
    dirichletEnergy=True,
)

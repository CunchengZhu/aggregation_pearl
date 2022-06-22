import pymem3dg.visual as dg_vis
import pymem3dg.read as dg_read
import matplotlib.pyplot as plt
import numpy as np
import polyscope as ps
import pymem3dg as dg
import sys
####################################################
#                 Initialize pathes                #
####################################################
folder = "../.." # run8
var = 80
subFolder = folder + f"/chi{var}"

sys.path.append(folder)
import parameters as para
xi, A_bar, R_bar, Kb, h = para.scalingVariables()
parameters = para.parameters(xi, A_bar, R_bar, Kb)
parameters.aggregation.chi = var * Kb / R_bar**2
trajNc = subFolder + "/traj.nc"
dg_vis.animate(
    trajNc,
    frames=np.arange(0, dg_read.sizeOf(trajNc) - 1, 1),
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
    entropyPotential=True,
    inPlaneFluxForm=True,
)
ps.show()
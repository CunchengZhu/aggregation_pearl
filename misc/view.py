import pymem3dg as dg
import pymem3dg.visual as dg_vis
import pymem3dg.util as dg_util
import pymem3dg.read as dg_read
import polyscope as ps
import imp
import matplotlib.pyplot as plt
import numpy as np

####################################################
#                 Initialize pathes                #
####################################################
# folder = "/beads_on_pearl"
# folder = "/intermediary"
# folder = r"H:/Shared drives/Rangamani Lab Drive/Cuncheng Zhu/manuscript_in_process/feng_pearling/repo/aggregation_pearl/chi20"
folder = ".."
# folder = "/no_pressure"
# folder = "/tube"

####################################################
#                  visualize .ply                  #
####################################################
# ply =  folder + "f15_t15148_.ply"
# face, vertex = dg.readMesh(ply)
# # print(dg.readData(ply))
# # print(dg.readData(ply, 'vertex'))
# # H = dg.readData(ply, 'vertex', 'mean_curvature')
# # Fb = dg.readData(ply, 'vertex', 'bending_force')

# ps.init()
# ps_mesh = ps.register_surface_mesh("membrane", vertex, face)
# # ps_mesh.add_scalar_quantity("mean_curvature", H, enabled=True)
# # ps_mesh.add_scalar_quantity("bending_force", Fb, enabled=True, vminmax=(-1e-5, 1e-5))
# ps.set_up_dir("z_up")
# ps.show()

####################################################
#            visualize animated .nc                #
####################################################
parameterFile = imp.load_source("module.name", folder + "/parameters.py")
xi, Ksg, Kb, R_bar, h = parameterFile.scalingVariables()
parameters = parameterFile.parameters(xi, R_bar, Kb, Ksg)
dg_vis.animate(
    folder + "/traj.nc",
    parameters=parameters,
    meanCurvature=True,
    gaussianCurvature=True,
    bendingForce=True,
    externalForce=True,
    mechanicalForce=True,
    capillaryForce=True,
    lineCapillaryForce=True,
    osmoticForce=True,
    chemicalPotential=True,
    bendingPotential=True,
    diffusionPotential=True,
    aggregationPotential=True,
    adsorptionPotential=True,
    springForce=True,
    entropyForce=True)

####################################################
#            visualize energy                      #
####################################################
dg_vis.getEnergyTrajectoryFromNc(
    trajFile=folder + "/traj.nc",
    parameters=parameters,
    # potentialEnergy=True,
    # kineticEnergy=True,
    # totalEnergy=True,
    bendingEnergy=True,
    # externalWork=True,
    # deviatoricEnergy=True,
    # surfaceEnergy=True,
    pressureEnergy=True,
    # adsorptionEnergy=True,
    aggregationEnergy=True,
    entropyEnergy=True,
    # edgeSpringEnergy=True,
    # faceSpringEnergy=True,
    # lcrSpringEnergy=True,
    dirichletEnergy=True
)

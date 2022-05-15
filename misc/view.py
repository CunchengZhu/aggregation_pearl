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
# trajFile = folder + "/traj.nc"

# parameterFile = imp.load_source("module.name", folder + "/parameters.py")
# xi, R_bar, Kb, h = parameterFile.scalingVariables()
# p = parameterFile.parameters(xi, R_bar, Kb)

# frameLim = (0, dg_read.sizeOf(trajFile))
# frameNum = frameLim[1] - frameLim[0]
# time = np.zeros(frameNum)
# kineticEnergy = np.zeros(frameNum)
# potentialEnergy = np.zeros(frameNum)
# externalWork = np.zeros(frameNum)
# totalEnergy = np.zeros(frameNum)
# volume = np.zeros(frameNum)
# for frame in range(frameNum):
#     system = dg.System(
#         trajFile,
#         frame,
#         p
#     )
#     system.initialize(nMutation = 0, ifMute = True)
#     time[frame] = system.time
#     volume[frame] = system.volume
#     system.computeTotalEnergy()
#     energy = system.getEnergy()
#     kineticEnergy[frame] = energy.kineticEnergy
#     potentialEnergy[frame] = energy.potentialEnergy
#     if frame != 0:
#         externalWork[frame] = externalWork[
#             frame - 1
#         ] + system.computeIntegratedPower(time[frame] - time[frame - 1])
# totalEnergy = potentialEnergy + kineticEnergy - externalWork
# reducedVolume = volume / (3.14 * 4 / 3)

# # plotting
# fig, ax1 = plt.subplots()
# color = 'tab:red'
# ax1.set_xlabel('time')
# ax1.set_ylabel('energy', color=color)
# ax1.plot(time, totalEnergy, color=color)
# ax1.tick_params(axis='y', labelcolor=color)

# ax2 = ax1.twinx()
# color = 'tab:blue'
# ax2.set_ylabel('reduced volume', color=color)
# ax2.plot(time, reducedVolume, color=color)
# ax2.tick_params(axis='y', labelcolor=color)

# fig.tight_layout()
# plt.show()
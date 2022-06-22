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
# folder = r"H:/Shared drives/Rangamani Lab Drive/Cuncheng Zhu/manuscript_in_process/feng_pearling/repo/aggregation_pearl/chi20"
folder = "../run8"
var = 100
subFolder = folder + f"/chi{var}"
sys.path.append(folder)
import parameters as ps
xi, A_bar, R_bar, Kb, h = ps.scalingVariables()
parameters = ps.parameters(xi, A_bar, R_bar, Kb)
parameters.aggregation.chi = var * Kb / R_bar**2
trajNc = subFolder + "/traj.nc"
####################################################
#                  visualize .ply                  #
####################################################
# dg_vis.visualizePly(plyFile, meanCurvature=True, gaussianCurvature=True)

# ps.init()
# face, vertex = dg.getIcosphere(1, 4)
# dg_vis.polyscopeStyle(True)
# ps_mesh = ps.register_surface_mesh("mesh", vertex, face, smooth_shade=True)
# ps.show()

# dg_vis.animate(
#     trajNc,
#     frames=[dg_read.sizeOf(trajNc) - 1],
# )
# ps.screenshot("hello.png")
####################################################
#            visualize animated .nc                #
####################################################
# dg_vis.animate(
#     trajNc,
#     frames=[5],
#     parameters=parameters,
#     meanCurvature=True,
#     gaussianCurvature=True,
#     bendingForce=True,
#     externalForce=True,
#     mechanicalForce=True,
#     capillaryForce=True,
#     lineCapillaryForce=True,
#     osmoticForce=True,
#     springForce=True,
#     entropyForce=True,
#     chemicalPotential=True,
#     bendingPotential=True,
#     diffusionPotential=True,
#     aggregationPotential=True,
#     adsorptionPotential=True,
#     entropyPotential=True,
#     inPlaneFluxForm=True,
# )
# ps.show()
####################################################
#                         plots                    #
####################################################
sp_size = (3,2)
fig, axs = plt.subplots(sp_size[0], sp_size[1])
dg_vis.matplotlibStyle(9, 10, 12)
fig.set_size_inches(8, 9)
count = 0
dg_vis.plotEnergy(
    axs[np.unravel_index(count, sp_size, "F")],
    trajNc,
    parameters,
    # logScale=True,
    zeroing=True,
    potentialEnergy=True,
    # kineticEnergy=True,
    # totalEnergy=True,
    bendingEnergy=True,
    # externalWork=True,
    # deviatoricEnergy=True,
    # surfaceEnergy=True,
    # pressureEnergy=True,
    # adsorptionEnergy=True,
    aggregationEnergy=True,
    entropyEnergy=True,
    # edgeSpringEnergy=True,
    # faceSpringEnergy=True,
    # lcrSpringEnergy=True,
    dirichletEnergy=True,
)
count = count + 1

dg_vis.plotChemicalPotentials(
    axs[np.unravel_index(count, sp_size, "F")],
    trajNc,
    parameters,
    # logScale=True,
    bendingPotential=True,
    # deviatoricPotential=True,
    aggregationPotential=True,
    entropyPotential=True,
    # dirichletPotential=True,
    # adsorptionPotential=True,
)
count = count + 1

dg_vis.plotMechanicalForces(
    axs[np.unravel_index(count, sp_size, "F")],
    trajNc,
    parameters,
    # logScale=True,
    bendingForce=True,
    # capillaryForce=True,
    # osmoticForce=True,
    # adsorptionForce=True,
    lineCapillaryForce=True,
    aggregationForce=True,
    # externalForce=True,
    entropyForce=True,
    # springForce=True,
)
count = count + 1

# count = np.prod(sp_size) - 1
# dg_vis.overlayColorMap(
#     axs[np.unravel_index(count, sp_size, "F")],
#     "hello.png",
#     [0, 1],
#     "$\phi$",
#     orientation="vertical",
# )

fig.tight_layout()
# plt.show()
plt.savefig("hello.png")

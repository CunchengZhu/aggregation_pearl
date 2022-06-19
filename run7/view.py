import pymem3dg.visual as dg_vis
import matplotlib.pyplot as plt
from tqdm.contrib.concurrent import process_map
import parameters as ps
import numpy as np

def plotTrajectory(trajNc, parameters, figName):
    sp_size = (2, 2)
    fig, axs = plt.subplots(sp_size[0], sp_size[1])
    dg_vis.matplotlibStyle(9, 10, 12)
    fig.set_size_inches(8, 6)
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

    dg_vis.plotProteinDensity(
        axs[np.unravel_index(count, sp_size, "F")], trajNc, parameters, logScale=True
    )
    count = count + 1

    dg_vis.plotChemicalPotentials(
        axs[np.unravel_index(count, sp_size, "F")],
        trajNc,
        parameters,
        logScale=True,
        bendingPotential=True,
        # deviatoricPotential=True,
        aggregationPotential=True,
        entropyPotential=True,
        dirichletPotential=True,
        # adsorptionPotential=True,
    )
    count = count + 1

    dg_vis.plotMechanicalForces(
        axs[np.unravel_index(count, sp_size, "F")],
        trajNc,
        parameters,
        logScale=True,
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

    fig.tight_layout()
    # plt.show()
    plt.savefig(figName)

def worker(var):
    dir = f'chi{var}'
    xi, A_bar, R_bar, Kb, h = ps.scalingVariables()
    parameters = ps.parameters(xi, A_bar, R_bar, Kb)
    parameters.aggregation.chi = var * Kb / R_bar**2
    plotTrajectory(dir+"/traj.nc", parameters, dir+".png")

def runPlots():
    jobs = []
    for v in np.arange(0, 120, 20):
        jobs.append(v)
    process_map(worker, jobs, max_workers=12)

if __name__ == "__main__":
    runPlots()
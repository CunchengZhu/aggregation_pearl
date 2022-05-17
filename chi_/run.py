import os
from tqdm.contrib.concurrent import process_map
from driver import localRun
import parameters as ps
import numpy as np

def worker(args):
    dir = args[0]
    var = args[1]
    if not os.path.exists(dir):
        os.mkdir(dir)
    os.chdir(dir)
    xi, A_bar, R_bar, Kb, h = ps.scalingVariables()
    parameters = ps.parameters(xi, A_bar, R_bar, Kb)
    parameters.aggregation.chi = var * Kb / R_bar**2
    localRun(parameters, CONTINUE = True)
    

# for i, kb in enumerate(np.arange(0,0.1,0.01)):
#     for replicate in np.arange(0,10000):
#         path = f'run_{i}_{replicate}/'

def runSims():
    jobs = []
    for v in np.arrange(0, 10, 2):
        path = f'chi{v}'
        jobs.append((path, v))
    process_map(worker, jobs, max_workers=12)

if __name__ == "__main__":
    runSims()
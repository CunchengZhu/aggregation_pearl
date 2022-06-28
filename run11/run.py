import os
from tqdm.contrib.concurrent import process_map
from driver import localRun
import parameters as ps
import numpy as np

def worker(args):
    dir = args[0]
    var = args[1]
    CONTINUE = args[2]
    if not os.path.exists(dir):
        os.mkdir(dir)
    os.chdir(dir)
    xi, A_bar, R_bar, Kb, h = ps.scalingVariables()
    parameters = ps.parameters(xi, A_bar, R_bar, Kb)
    parameters.aggregation.chi = var * Kb / R_bar**2
    localRun(parameters, CONTINUE = CONTINUE)

def runSims(CONTINUE = True):
    jobs = []
    for v in np.arange(0, 80, 10):
        path = f'chi{v}'
        jobs.append((path, v, CONTINUE))
    process_map(worker, jobs, max_workers=12)

if __name__ == "__main__":
    runSims(CONTINUE=True)

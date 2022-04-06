import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import qmc
import scipy


def halton_generate_SND(n_time_steps, n_paths):
    halton_SND=np.zeros((n_time_steps,n_paths),np.float64)
    
    sampler = qmc.Halton(d=n_time_steps, scramble=True)
    Halton = sampler.random(n_paths)
    rand = scipy.stats.norm.ppf(Halton)
    
    for t in range(n_time_steps):
        halton_SND[t] = rand[:,t]
    return halton_SND
"""This script make fof output halo readble by covo to study correlations
Oriane Laurens
June 2025"""

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import pandas as pd
from astropy import units as u
import math
mp.rcParams['agg.path.chunksize'] = 10000 # to adjust size (suggested by a previous error)

# READING THE COVO OUTPUT FILE
file_name = f"/path/to/file" #koki
print("About to read {}".format(file_name))
data = pd.read_csv(file_name, sep='\s+')

# GET CORRELATION DATA
def plot_data_correl():
    shift = .5
    r = data.iloc[:, 1] #distance
    r12_v1a = data.iloc[:, 3] - shift #correlation vector a shifted
    r12_v1a_std = data.iloc[:, 4] #standard deviation correlation vector a
    r12_v1b= data.iloc[:, 5] - shift#correlation vector b
    r12_v1b_std = data.iloc[:, 6] #standard deviation correlation vector b
    fig, ax = plt.subplots(2, 1, figsize=(8,5))
    ax[0].errorbar(r, r12_v1a, yerr = r12_v1a_std, xerr = None, label = 'vector a', color='blue')
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Distance r')
    ax[0].set_ylabel('Correlation')
    ax[0].set_title('Correlation function of semi-major axis a (run394, z=0)')
    ax[0].legend()
    ax[0].grid(True)
    ax[1].errorbar(r, r12_v1b, yerr = r12_v1b_std, xerr = None, label = 'moment of inertia', color='blue')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Distance r')
    ax[1].set_ylabel('Correlation')
    ax[1].set_title('Correlation function of moment of inertia (run394, z=0)')
    ax[1].legend()
    ax[1].grid(True)
    plt.tight_layout()
    plt.savefig(f"correlation_functions_-.5_z0_run394_plot.png", dpi = 300)

# EXECUTION
plot_data_correl()


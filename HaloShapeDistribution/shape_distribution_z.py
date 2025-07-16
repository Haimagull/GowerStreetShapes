"""
This script aims to evaluate DM halo shape distribution dependency on redshift variation. 
There is a limit number of particles of a halo set to 300, we only keep those with Np_halo > 300,
based on the suggestions of https://doi.org/10.48550/arXiv.1203.6833 .
One of the main goals is to study correlations between shapes.

Oriane Laurens

June 2025
"""

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from astropy import constants as const
from astropy import units as u
import math
mp.rcParams['agg.path.chunksize'] = 10000 # to adjust size (suggested by a previous error)

# CALCULATING PARTICLE MASS
def run_control_par(run) :# reading control.par file from run
    #file_name = f"/Users/orianenyembo/Internship_2025_UCL/control.par"
    file_name = f"/share/testde/ucapwhi/GowerStreetExtendedSims/runsW/{run}/control.par"
    variable = {}
    print("About to read {}".format(file_name))
    with open(file_name, "r") as in_file :
        exec(in_file.read(), {"math": math}, variable)
    dBoxSize = variable["dBoxSize"]
    nGrid = variable["nGrid"]
    dOmega0 = variable["dOmega0"]
    h = variable["h"]
    return dBoxSize, nGrid, dOmega0, h


# READING THE FOF FILE

def fof_file_format_experiment(time_slice) :
    fof_type = np.dtype([
       ('position_of_deepest_potential', np.float32, (3,)),
       ('deepest_potential', np.float32, 1),
       ('shrinking_sphere_centre', np.float32, (3,)),
       ('position_of_centre_of_mass', np.float32, (3,)),
       ('velocity_of_centre_of_mass', np.float32, (3,)),
       ('angular_momentum', np.float32, (3,)),
       ('moment_of_inertia', np.float32, (6,)),
       ('velocity_dispersion', np.float32, 1),
       ('r_max', np.float32, 1),
       ('mass', np.float32, 1),
       ('mass_env_0', np.float32, 1),
       ('mass_env_1', np.float32, 1), 
       ('half_mass_radius', np.float32, 1)])

    #file_name = f"/export/koki1/visitor1/Gowersims/runs/runsW/run243fof/fof/run.{time_slice}.fofstats.0"
    #file_name = f"/Users/orianenyembo/Internship_2025_UCL/GowerStreetExtendedSims/runsW/run243/fof/run.{time_slice}.fofstats.0"
    file_name = f"/share/testde/ucapwhi/GowerStreetExtendedSims/runsW/{run}/fof/run.{time_slice}.fofstats.0"
    print("About to read {}".format(file_name))
    with open(file_name, "rb") as in_file:
        data = np.fromfile(in_file, dtype = fof_type)
    print(data['position_of_deepest_potential'])
    print(data['position_of_centre_of_mass'])
    print(data.shape)
    return data

def get_columns(time_slice) :
    data = fof_file_format_experiment(time_slice)
    position_dp = data['position_of_deepest_potential'] # getting columns
    dp = data['deepest_potential']
    ssc = data['shrinking_sphere_centre']
    position_cdm = data['position_of_centre_of_mass']
    velocity_cdm = data['velocity_of_centre_of_mass']
    L_ang = data['angular_momentum']
    inertia = data['moment_of_inertia'] #'loading' inertia tensor file data
    sigma_disp = data['velocity_dispersion']
    r_max = data['r_max']
    mass = data['mass']
    mass_env_0 = data['mass_env_0']
    mass_env_1 = data['mass_env_1']
    hmr = data['half_mass_radius']
    return position_dp, dp, ssc, position_cdm, velocity_cdm, L_ang, inertia, sigma_disp, r_max, mass, mass_env_0, mass_env_1, hmr

def particle_mass(time_slice, run) :
    """This function gives particle mass"""
    # Get halo mass
    data = fof_file_format_experiment(time_slice) 
    mass = data['mass']
    # Get parameters from control.par (dBoxSize in Mpc/h)
    dBoxSize, nGrid, dOmega0, h = run_control_par(run)
    # Calc.
    G_c = 1#const.G #m3 / (kg s2); Gravi.const. = 6.6743e-11 w/ uncertainty = 1.5e-15
    H = h*100*(u.km / u.s / u.Mpc)
    rho_c = (3*H**2)/(8*np.pi*(G_c*(u.Mpc**3 / (u.kg * u.s**2))))#.to(u.Mpc**3 / (u.kg * u.s**2)))) #100 units : km.s-1.Mpc-1
    M_box = (dOmega0 * rho_c * (dBoxSize*(u.Mpc/h))**3).to(u.kg) #Box mass in kg (to kg cause units km2/Mpc2)
    M_part = (M_box/nGrid**3).to('Msun') #Particle mass in solar mass Msun converted from kg w/ astropy units
    return M_part #Particle mass in Msun

# GET RATIOS FROM SEMI-AXIS AND INERTIA MATRIX
def axis_ratio(time_slice, run) :
    data = fof_file_format_experiment(time_slice) 
    mass = data['mass']*u.Msun
    inertia = data['moment_of_inertia']
    M_part = particle_mass(time_slice, run)
    r_ca_data=[] #list of c/a values
    r_ba_data=[] #list of b/a values
    for i,row in enumerate(inertia) :
        Np_halo = mass[i] / M_part #Nb of part. in halo
        if Np_halo > 300 :
            """Construction of a, b, c from inertia matrix
            We have 6 independent components
            Diagonal terms : d1=Ixx; d2=Iyy; d3=Izz
            Triangle terms : t1=Ixy=Iyx; t2=Ixz=Izx; t3=Iyz=Izy
            + Output is two ratios, first c/a then b/a"""
            d1, t1, t2, d2, t3, d3 = row #get matrix independent components from file
            I = np.array([[d1,t1,t2], [t1,d2,t3], [t2,t3,d3]]) #reconstructing inertia matrix
            eigenval, eigenvect = np.linalg.eigh(I) #get eigen values sorted with a>=b>=c
            c, b, a = np.sqrt(abs(eigenval))  # 
            r_ca_row = c/a
            r_ba_row = b/a
            r_ca_data.append(r_ca_row)
            r_ba_data.append(r_ba_row)
        else :
            continue
    return r_ca_data, r_ba_data

# PLOT
def plot_contour(time_slice, run) :
    r_ca_data, r_ba_data = axis_ratio(time_slice, run)
    xmin, xmax, ymin, ymax = 0, 1, 0, 1 #(normalized by a)
    hist, xedges, yedges = np.histogram2d(r_ca_data, r_ba_data, bins=100, range=[[xmin, xmax], [ymin, ymax]])
    xcenters = (xedges[:-1] + xedges[1:])/2 #centered bin
    ycenters = (yedges[:-1] + yedges[1:])/2
    X, Y = np.meshgrid(xcenters, ycenters)
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = hist.T #"density"
    #plot
    plt.figure(figsize=(8, 6))
    plt.title(f'Distribution of halo shape in {run} ({time_slice}/00100)')
    plt.xlabel('c/a')
    plt.ylabel('b/a')
    plt.contour(X, Y, Z, levels = 5, origin='lower')
    plt.colorbar(label='Density')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.tight_layout()
    plt.savefig(f"hsd_{run}.{time_slice}_contour.png")
    plt.close()
    r_ca_data.clear()
    r_ba_data.clear()

# EXECUTION

runs = ['run394','run273','run378','run355']
#times = [str(i).zfill(5) for i in range(1, 101)]
#for time_slice in times :
for run in runs :
    plot_contour(time_slice='00100', run=run)

""" If you want to run it on several runs add :
#runsW = [str(i).zfill(3) for i in range(1, 897)]
#for run in runsW :
#    times = [str(i).zfill(5) for i in range(1, 101)]
#    for time_slice in times :
#        plot_contour(time_slice)
#in the execution part instead of run = "run243" """


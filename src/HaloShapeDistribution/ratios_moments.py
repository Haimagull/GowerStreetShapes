"""
This calculates the moments (mean, variance, skewness) over redshifts for each run. 
If you already calculated the shape distribution for each time slice of a run, this is redundant.
It saves first order moments of axis ratios in a csv file with the assiocated redshift.
The output csv file can then be used to calculate error bars and plot mean, variance and skewness 
as a function of redshift for chosen runs.

Warning : run ??? must be called with 'run???' 

Oriane Laurens

July 2025
"""

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from astropy import constants as const
from astropy import units as u
from scipy.stats import skew
import math
mp.rcParams['agg.path.chunksize'] = 10000 # to adjust size (suggested by a previous error)

import csv
import os

# CALCULATING PARTICLE MASS
def run_control_par(control_path):#run) :# reading control.par file from run
    #file_name = f"/share/testde/ucapwhi/GowerStreetExtendedSims/runsW/{run}/control.par"
    file_name = control_path
    variable = {}
    print("About to read {}".format(file_name))
    with open(file_name, "r") as in_file :
        exec(in_file.read(), {"math": math}, variable)
    dBoxSize = variable["dBoxSize"]
    nGrid = variable["nGrid"]
    dOmega0 = variable["dOmega0"]
    h = variable["h"]
    return dBoxSize, nGrid, dOmega0, h

# CALCULATION 1ST TO 3RD MOMENTS (MEAN, VARIANCE, SKEWNESS)
def moments(ratio) : #moments(time_slice, run, ratio) : #ratio is r_ca_data or r_ba_data from axis_ratio(time_slice, run)
    mean = np.mean(ratio) # mean, 1st moment
    variance = np.var(ratio, ddof=1) # variance, 2nd moment
    skewness = skew(ratio) # skewness, 3rd moment
    return mean, variance, skewness
    
# READING REDSHIFT VALUES
def read_z_values(z_values_path):#run):
    file_name = z_values_path #f'/share/testde/ucapwhi/GowerStreetExtendedSims/runsW/{run}/z_values.txt'
    data = np.genfromtxt(file_name, delimiter=',')#, names=True, dtype=None, encoding='utf-8')
    return data[:, 2]
  
# READING THE FOF FILE
def fof_file_format_experiment(time_slice, run, runs_path) :
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
    sub_path = f"/{run}/fof/run.{time_slice}.fofstats.0"
    file_name = runs_path + sub_path #f"/share/testde/ucapwhi/GowerStreetExtendedSims/runsW/{run}/fof/run.{time_slice}.fofstats.0"
    print("About to read {}".format(file_name))
    with open(file_name, "rb") as in_file:
        data = np.fromfile(in_file, dtype = fof_type)
    print(data.shape)
    return data

def particle_mass(run, data, control_path) :
    """This function gives particle mass"""
    # Get halo mass
    mass = data['mass']
    # Get parameters from control.par (dBoxSize in Mpc/h)
    dBoxSize, nGrid, dOmega0, h = run_control_par(control_path)
    # Calc.
    G_c = 1#const.G #m3 / (kg s2); Gravi.const. = 6.6743e-11 w/ uncertainty = 1.5e-15
    H = h*100*(u.km / u.s / u.Mpc)
    rho_c = (3*H**2)/(8*np.pi*(G_c*(u.Mpc**3 / (u.kg * u.s**2))))#.to(u.Mpc**3 / (u.kg * u.s**2)))) #100 units : km.s-1.Mpc-1
    M_box = (dOmega0 * rho_c * (dBoxSize*(u.Mpc/h))**3).to(u.kg) #Box mass in kg (to kg cause units km2/Mpc2)
    M_part = (M_box/nGrid**3).to('Msun') #Particle mass in solar mass Msun converted from kg w/ astropy units
    return M_part #Particle mass in Msun

# GET RATIOS FROM SEMI-AXIS AND INERTIA MATRIX
def axis_ratio(run, data, control_path) :
    mass = data['mass']*u.Msun
    inertia = data['moment_of_inertia']
    M_part = particle_mass(run, data, control_path)
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
            a_2, b_2, c_2 = np.sort(eigenval)[::-1] #decreasing sorting so a, b, c
            a = np.sqrt(np.abs(a_2))
            b = np.sqrt(np.abs(b_2))
            c = np.sqrt(np.abs(c_2))
            #c, b, a = np.sqrt(abs(eigenval))  # 
            r_ca_row = c/a
            r_ba_row = b/a
            r_ca_data.append(r_ca_row)
            r_ba_data.append(r_ba_row)
        else :
            continue
    return r_ca_data, r_ba_data

# PLOT MOMENTS
def plot_moments(times, run, z_values_path) :
    redshift_val = read_z_values(z_values_path)
    if len(times) != len(redshift_val):
        raise ValueError(f"{len(times)} time slices but {len(redshift_val)} redshifts, must be the same length")
    means_r_ca, vars_r_ca, skews_r_ca = [], [], []
    means_r_ba, vars_r_ba, skews_r_ba = [], [], []
    means_r_ca_err, vars_r_ca_err, skews_r_ca_err = [], [], []
    means_r_ba_err, vars_r_ba_err, skews_r_ba_err = [], [], []
    #output directory to save results
    output_dir = os.path.join("results", run, "moments")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"moments_{run}.csv")
    #save everything in a csv file
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["TimeSlice", "Redshift", "Mean_c/a", "Var_c/a", "Skew_c/a", "Mean_b/a", "Var_b/a", "Skew_b/a"])
        for time_slice, z in zip(times, redshift_val) :
            data = fof_file_format_experiment(time_slice, run, runs_path) 
            r_ca_data, r_ba_data = axis_ratio(run, data, control_path)
            mean_r_ca, var_r_ca, skew_r_ca = moments(r_ca_data)
            mean_r_ba, var_r_ba, skew_r_ba = moments(r_ba_data)
            means_r_ca.append(mean_r_ca)
            vars_r_ca.append(var_r_ca)
            skews_r_ca.append(skew_r_ca)
            means_r_ba.append(mean_r_ba)
            vars_r_ba.append(var_r_ba)
            skews_r_ba.append(skew_r_ba)
            #write to csv
            writer.writerow([time_slice, z, mean_r_ca, var_r_ca, skew_r_ca, mean_r_ba, var_r_ba, skew_r_ba])

# EXECUTION
print("Warning : Simulations fof outputs of each run must be organised the following way : /run/fof/run.time_slice.fofstats.0")
#print("If you")
control_path = str(input("Please, indicate the path to your control.par file :"))
z_values_path = str(input(("Please, indicate the path to your z_values.txt file :")))
runs_path = str(input("Please, indicate the path to your runs :"))
print("Please, indicate the names of your chosen runs (e.g. [run243, run344, run566]), when you finish type 'stop' and enter.")
runs = []
var = ""
while var != "stop" :
    var = input("->")
    if var == "" :
        continue
    runs.append(var)
runs.pop()
    
for run in runs :
    times = [str(i).zfill(5) for i in range(1, 101)] #loops inside the function
    plot_moments(times, run, z_values_path)


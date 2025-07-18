"""
This script calculates semi axis a, b, c of halos (min. 300 particles) 
found in the fof output of Gower Street Simulations.
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
import math
mp.rcParams['agg.path.chunksize'] = 10000 # to adjust size (suggested by a previous error)
import os
import csv

# CALCULATING PARTICLE MASS
def run_control_par(run) :# reading control.par file from run
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

def read_z_values(run):
    file_name = f'/share/testde/ucapwhi/GowerStreetExtendedSims/runsW/{run}/z_values.txt'
    data = np.genfromtxt(file_name, delimiter=',')#, names=True, dtype=None, encoding='utf-8')
    return data[:, 2]

# READING THE FOF FILE

def fof_file_format_experiment(time_slice, run) :
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
    file_name = f"/share/testde/ucapwhi/GowerStreetExtendedSims/runsW/{run}/fof/run.{time_slice}.fofstats.0"
    print("About to read {}".format(file_name))
    with open(file_name, "rb") as in_file:
        data = np.fromfile(in_file, dtype = fof_type)
    print(data['position_of_deepest_potential'])
    print(data['position_of_centre_of_mass'])
    print(data.shape)
    return data

def particle_mass(run, data) :
    """This function gives particle mass"""
    # Get halo mass
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
def shape(run, data, time_slice, redshift_val) : 
    mass = data['mass']*u.Msun
    inertia = data['moment_of_inertia']
    M_part = particle_mass(run, data)
    a_data = [] #list of a values
    b_data = [] #list of b values
    c_data = [] #list of c values
    eigenvect_data = []
    output_dir = os.path.join("results", run, "moments") #output directory to save results
    os.makedirs(output_dir, exist_ok=True) #check it exists
    output_file = os.path.join(output_dir, f"shape_{run}.csv") #output file
    with open(output_file, mode='w', newline='') as file : #save everything in a csv file
        writer = csv.writer(file)
        writer.writerow(["Time Slice", "Redshift", "a", "b", "c"])
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
            a_data.append(a)
            b_data.append(b)
            c_data.append(c)
            #write to csv
            writer.writerow([time_slice, redshift_val, a_data, b_data, c_data])
        else :
            continue
    return a_data, b_data, c_data

def get_run_shapes(run) :
    times = [str(i).zfill(5) for i in range(1, 101)]
    redshift_val = read_z_values(run)
    for time_slice,z in zip(times, redshift_val) :
        data = fof_file_format_experiment(time_slice, data)
        shape(run, data, time_slice=time_slice, redshift_val=redshift_val)

# EXECUTION
runs = ['run???', 'run???'] #replace with wanted runs
for run in runs :
    get_run_shapes(run=run)



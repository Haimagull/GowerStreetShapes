"""
This script make fof output halo readble by covo to study correlations. 
In this version, we still have a minimum of 300 particles for a halo and a
selection has also been made on halos mass to have a smaller catalog to
process, we only keep the 10% most massive halos.

Oriane Laurens

June 2025
"""

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import pandas as pd
from astropy import units as u
import math
mp.rcParams['agg.path.chunksize'] = 10000 # to adjust size (suggested by a previous error)

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
    file_name = runs_path + sub_path
    print("About to read {}".format(file_name))
    with open(file_name, "rb") as in_file:
        data = np.fromfile(in_file, dtype = fof_type)
    print(data.shape)
    return data

# GET BOX SIZE
def box_size_function(run, control_path) :# reading control.par file from run
    file_name = control_path
    variable = {}
    print("About to read {}".format(file_name))
    with open(file_name, "r") as in_file :
        exec(in_file.read(), {"math": math}, variable)
    dBoxSize = variable["dBoxSize"]
    return dBoxSize

# GET RATIOS FROM SEMI-AXIS AND INERTIA MATRIX
def axis_vect(time_slice, run, runs_path) :
    data = fof_file_format_experiment(time_slice, run, runs_path)
    inertia = data['moment_of_inertia'] #'loading' inertia tensor file data
    eigenvects = []
    mass = data['mass']
    quantile = np.percentile(mass, 90) #90th percentile
    for mass_row, inertia_row in zip(mass, inertia) :
        """Construction of a, b, c from inertia matrix
        We have 6 independent components
        Diagonal terms : d1=Ixx; d2=Iyy; d3=Izz
        Triangle terms : t1=Ixy=Iyx; t2=Ixz=Izx; t3=Iyz=Izy
        + Output is two ratios, first c/a then b/a"""
        if mass_row > quantile : #only selects above the 90-th centile so only the top 10%
            d1, t1, t2, d2, t3, d3 = inertia_row #get matrix independent components from file
            I = np.array([[d1,t1,t2], [t1,d2,t3], [t2,t3,d3]]) #reconstructing inertia matrix
            eigenval, eigenvect = np.linalg.eigh(I) #get eigen values sorted with a>=b>=c
            eigenvects.append(eigenvect)
    return eigenvects

# STORE DATA FOR COVO
def get_data(time_slice, run, runs_path, control_path):
    data = fof_file_format_experiment(time_slice, run, runs_path)
    position_cdm = data['position_of_centre_of_mass']
    inertia_vect = data['moment_of_inertia']
    box_size = box_size_function(run, control_path)
    eigenvects = axis_vect(time_slice, run, runs_path)
    xs = []
    ys = []
    zs = []
    axs = []
    ays = []
    azs = []
    bxs = []
    bys = []
    bzs = []
    for i, row in enumerate(eigenvects):
        xs.append(position_cdm[i, 0]*box_size)
        ys.append(position_cdm[i, 1]*box_size)
        zs.append(position_cdm[i, 2]*box_size)
        eigenvect_a = row[0,:] #get first vector, here it is semi-minor axis c. If you want b use [1,:] or a use [2,:]
        inertia = inertia_vect[i, :]
        ax = eigenvect_a[0]
        axs.append(ax)
        ay = eigenvect_a[1]
        ays.append(ay)
        az = eigenvect_a[2]
        azs.append(az)
        bx = inertia[0]
        bxs.append(bx)
        by = inertia[1]
        bys.append(by)
        bz = inertia[2]
        bzs.append(bz)
    halo_data = {
      'x' : xs,
      'y' : ys,
      'z' : zs,
      'ax' : axs,
      'ay' : ays,
      'az' : azs,
      'bx' : bxs,
      'by' : bys,
      'bz' : bzs
   }
    data_file = pd.DataFrame(halo_data)
    data_file.to_csv(f'covo_typefile_{run}_{time_slice}_aI.csv', index=False, header=False)

# EXECUTION
print("Warning : Simulations fof outputs of each run must be organised the following way : /run/fof/run.time_slice.fofstats.0")
control_path = str(input("Please, indicate the path to your control.par file :"))
runs_path = str(input("Please, indicate the path to your runs :"))
print("Please, indicate the names of your chosen runs (e.g. [run243, run344, run566]), when you finish type 'stop' and enter.")
run = str(input("Please indicate the chosen run (e.g. run243) :"))
time_slice = str(input("Please indicate the chosen time slice :")).zfill(5)
get_data(time_slice=time_slice, run=run, runs_path=runs_path, control_path=control_path)

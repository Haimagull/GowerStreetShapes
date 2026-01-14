"""
This script make fof output halo readble by IAcorr to study 
correlations of ellipticity for the ten percent most massive halos
Oriane Laurens
December 2025
"""

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import pandas as pd
from astropy import units as u
from astropy.coordinates import CartesianRepresentation, SkyCoord
import math
mp.rcParams['agg.path.chunksize'] = 10000 # to adjust size (suggested by a previous error)

# PATHS
def z_values_file(run) :
 return '/Users/orianenyembo/Internship_2025_UCL/ellcorr/run243/z_values.txt'
def control_file(run) :
    return '/Users/orianenyembo/Internship_2025_UCL/ellcorr/run243/control.par'
def fof_run_file(time_slice) :
    return '/Users/orianenyembo/Internship_2025_UCL/ellcorr/run243/run.00056.fofstats.0'

# READING THE FOF FILE

def fof_file_format_experiment(time_slice, runs_path) :
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
    print(data['position_of_deepest_potential'])
    print('data shape is',data.shape)
    return data

# GET BOX SIZE
def box_size_function(control_path) :# reading control.par file from run
    file_name = control_path
    variable = {}
    print("About to read {}".format(file_name))
    with open(file_name, "r") as in_file :
        exec(in_file.read(), {"math": math}, variable)
    dBoxSize = variable["dBoxSize"]
    return dBoxSize

# GET REDSHIFT
def read_z_values(z_values_path):
    file_name = z_values_path
    data = np.genfromtxt(file_name, delimiter=',')#, names=True, dtype=None, encoding='utf-8')
    return data[:, 2]

def twoD_projection(a, b, c, x, y, z) : #los is a string choice of line of sight
    """
    This code aims to project 3D ellipsoids on a 2D plane
    """
    los = np.array([x, y, z]) / np.sqrt(x*x + y*y + z*z)
    a_val, b_val, c_val = np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)
    #we use latest halo with eigenvect to get s=transpose(sxmu,symu,sLOSmu) vector
    s3_los = np.dot((c/c_val), los) #c
    mu3 = (c/c_val) - s3_los * los
    s2_los = np.dot((b/b_val), los) #b
    mu2 = (b/b_val) - s2_los * los
    s1_los = np.dot((a/a_val), los) #a
    mu1 = (a/a_val) - s1_los * los
    s1 = mu1.reshape(-1,1)
    s2 = mu2.reshape(-1,1)
    s3 = mu3.reshape(-1,1)
    #print(f'Unit eigenvectors of halo inertia tensor are : s1={s1}, s2={s2}, s3={s3}.')
    #get perpendicular s, sp
    sp1 = np.delete(mu1,2).reshape(-1,1)
    sp2 = np.delete(mu2,2).reshape(-1,1)
    sp3 = np.delete(mu3,2).reshape(-1,1)
    #get its transpose, sp_t
    sp1_t = sp1.reshape(1,-1)
    sp2_t = sp2.reshape(1,-1)
    sp3_t = sp3.reshape(1,-1)
    #get alpha
    alpha = (s1_los/a_val)**2 + (s2_los/b_val)**2 + (s3_los/c_val)**2
    #get k and transposed k_t
    k = (s1_los*sp1)/a_val**2 +(s2_los*sp2)/b_val**2 + (s3_los*sp3)/c_val**2
    k_t = k.reshape(1,-1)
    #get W_inv
    W_inv = ((sp1*sp1_t/a_val**2)-(k*k_t/alpha**2)) + ((sp2*sp2_t/b_val**2)-(k*k_t/alpha**2)) + ((sp3*sp3_t/c_val**2)-(k*k_t/alpha**2))
    """
    Calculating 2D ellipticity e1 + i e2 = ec_mod * exp(theta)
    """
    e1 = (W_inv[0,0]-W_inv[1,1])/(W_inv[0,0]+W_inv[1,1])
    e2 = 2*W_inv[0,1]/(W_inv[0,0]+W_inv[1,1])
    ec = complex(e1,e2) #complex ellipticity
    #ec_mod = abs(ec) #module of ellipticity
    phi_e = np.angle(ec)/2 #phase phi = theta/2 (cause spin two shape) with theta the argument of complex number ec
    return e1, e2, phi_e

#FUNCTION THAT TRANSFORMS X,Y,Z TO RA/DEC (with astropy)
def box_to_radec(x, y, z):
    cartesian = CartesianRepresentation(x*u.Mpc, y*u.Mpc, z*u.Mpc) # define cartesian
    proj_sky = SkyCoord(cartesian, frame='icrs') # convert in ra/dec, project on a sky
    return proj_sky.ra.deg, proj_sky.dec.deg #ra/dec in degrees !!

# GET RATIOS FROM SEMI-AXIS AND INERTIA MATRIX
def get_data_axis_vect(time_slice,run, runs_path, control_path, z_values_path) :
    z_values = read_z_values(z_values_path)
    z_val = z_values[int(time_slice)]
    data = fof_file_format_experiment(time_slice, runs_path)
    inertia = data['moment_of_inertia'] #'loading' inertia tensor file data
    position_cdm = data['position_of_centre_of_mass']
    box_size = box_size_function(run, control_path)
    ras = []
    decs = []
    zs = []
    eps1 = []
    eps2 = []
    masses = []
    mass = data['mass']
    quantile = np.percentile(mass, 90) #90th percentile
    for i, (mass_row, inertia_row) in enumerate(zip(mass, inertia)) :
        """
        Construction of a, b, c from inertia matrix
        We have 6 independent components
        Diagonal terms : d1=Ixx; d2=Iyy; d3=Izz
        Triangle terms : t1=Ixy=Iyx; t2=Ixz=Izx; t3=Iyz=Izy
        + Output is two ratios, first c/a then b/a
        """
        if mass_row > quantile : #only selects above the 90-th centile so only the top 10%
            d1, t1, t2, d2, t3, d3 = inertia_row #get matrix independent components from file
            I = np.array([[d1,t1,t2], [t1,d2,t3], [t2,t3,d3]]) #reconstructing inertia matrix
            eigenval, eigenvect = np.linalg.eigh(I) #get eigen values sorted with a>=b>=c
            a_vect = np.sqrt(eigenval[2]) * eigenvect[:,2]
            b_vect = np.sqrt(eigenval[1]) * eigenvect[:,1]
            c_vect = np.sqrt(eigenval[0]) * eigenvect[:,0]
            x = position_cdm[i, 0]*box_size
            y = position_cdm[i, 1]*box_size
            z = position_cdm[i, 2]*box_size #coordinate
            e1, e2, phi_e = twoD_projection(a_vect, b_vect, c_vect, x, y, z)
            eps1.append(e1)
            eps2.append(e2)
            ra, dec = box_to_radec(x, y, z)
            ras.append(ra)
            decs.append(dec)
            zs.append(z_val) #redshift
            masses.append(mass_row)
    halo_data = {
      'ra' : ras,
      'dec' : decs,
      'z' : zs,
      'eps1' : eps1,
      'eps2' : eps2,
      'mass' : masses
   }
    data_file = pd.DataFrame(halo_data)
    data_file.to_csv(f'IACorr_typefile_{time_slice}.csv', index=False, header=False)

#EXECUTION
print("Warning : Simulations fof outputs of each run must be organised the following way : /run/fof/run.time_slice.fofstats.0")
control_path = str(input("Please, indicate the path to your control.par file :"))
z_values_path = str(input(("Please, indicate the path to your z_values.txt file :")))
runs_path = str(input("Please, indicate the path to your runs :"))
print("Please, indicate the names of your chosen runs (e.g. [run243, run344, run566]), when you finish type 'stop' and enter.")
run = str(input("Please indicate the chosen run (e.g. run243) :"))
time_slice = str(input("Please indicate the chosen time slice :")).zfill(5)
get_data_axis_vect(time_slice=time_slice, run=run, runs_path=runs_path, control_path=control_path, z_values_path=z_values_path)

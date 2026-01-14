"""
HALO SHAPER
This code aims to provide a tool to visualize 3D ellipsoids (modelizing DM haloes)
and the impact of its semi-axis values. Also, there is a projection code from ellipsoid 
to 2D ellipse visualised by its complex ellipticity.

Oriane Laurens 

June 2025
"""

#Imports

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 
plt.ion() #interactive mode
plt.switch_backend('Qt5Agg') #graphic display method
from matplotlib.widgets import Slider, Button, RadioButtons #for widgets

#Starting figure
fig=plt.figure(figsize=(8,6))
ax=fig.add_subplot(projection='3d')
ax.set_position([0, 0.15, 1, 0.8])
 
#Initial parameters
(a0,b0,c0)=(1, 1, 1) #initially we have a sphere

#Halo parametrization
def halo(a, b, c) :
    """
    Parametrization of a 3D halo from semi-axis
    """
    theta = np.linspace(0,np.pi)
    phi = np.linspace(0,2*np.pi)
    theta, phi = np.meshgrid(theta, phi)
    x = a*np.sin(theta)*np.cos(phi)
    y = b*np.sin(theta)*np.sin(phi)
    z = c*np.cos(theta)
    return x, y, z

#2D Projection
#Constructing inertia matrix from a, b, c
from scipy.stats import ortho_group #to generate a random orthogonal matrix
def inertia_matrix(a, b, c) :
   diag_matrix = np.array([[a**2,0,0],[0,b**2,0],[0,0,c**2]]) #diagonal matrix of inertia
   #rot_matrix = ortho_group.rvs(3) #rotation matrix (random orthogonal matrix) #can be used to add a random alignment and see its impact on projection
   #rot_trans = np.transpose(rot_matrix) #transposed rotation matrix
   #first_mult = np.matmul(rot_matrix, diag_matrix)
   inertia_matrix = diag_matrix #np.matmul(first_mult, rot_trans)
   return inertia_matrix

def twoD_projection(a, b, c, los) : #los is a string choice of line of sight
    """
    This code aims to project 3D ellipsoids on a 2D plane
    """
    eigenval, eigenvect = np.linalg.eigh(inertia_matrix(a, b, c)) #get eigen values sorted with c=<b=<a
    #we use latest halo with eigenvect to get s=transpose(sxmu,symu,sLOSmu) vector
    if los == 'x' :
        mu3 = np.delete(eigenvect[0,:],2) #c
        s3_los = eigenvect[0,:][0]
        mu2 = np.delete(eigenvect[1,:],2) #b
        s2_los = eigenvect[1,:][0]
        mu1 = np.delete(eigenvect[2,:],2) #a
        s1_los = eigenvect[2,:][0]
        mu1 = np.append(mu1,s1_los)
        mu2 = np.append(mu2,s2_los)
        mu3 = np.append(mu3,s3_los)
    if los == 'y' :
        mu3 = np.delete(eigenvect[0,:],2) #c
        s3_los = eigenvect[0,:][1]
        mu2 = np.delete(eigenvect[1,:],2) #b
        s2_los = eigenvect[1,:][1]
        mu1 = np.delete(eigenvect[2,:],2) #a
        s1_los = eigenvect[2,:][1]
        mu1 = np.append(mu1,s1_los)
        mu2 = np.append(mu2,s2_los)
        mu3 = np.append(mu3,s3_los)
    if los == 'z' :
        mu3 = eigenvect[0,:] #vector for component 3
        mu2 = eigenvect[1,:] #idem 2
        mu1 = eigenvect[2,:] #idem 1
        s3_los = eigenvect[0,:][2]
        s2_los = eigenvect[1,:][2]
        s1_los = eigenvect[2,:][2]
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
    alpha = (s1_los/a)**2 + (s2_los/b)**2 + (s3_los/c)**2
    #get k and transposed k_t
    k = (s1_los*sp1)/a**2 +(s2_los*sp2)/b**2 + (s3_los*sp3)/c**2
    k_t = k.reshape(1,-1)
    #get W_inv
    W_inv = ((sp1*sp1_t/a**2)-(k*k_t/alpha**2)) + ((sp2*sp2_t/b**2)-(k*k_t/alpha**2)) + ((sp3*sp3_t/c**2)-(k*k_t/alpha**2))
    """Calculating 2D ellipticity"""
    e1 = (W_inv[0,0]-W_inv[1,1])/(W_inv[0,0]+W_inv[1,1])
    e2 = 2*W_inv[0,1]/(W_inv[0,0]+W_inv[1,1])
    ec = complex(e1,e2) #complex ellipticity
    ec_mod = abs(ec) #module of ellipticity
    phi_e = np.angle(ec)/2 #phase phi = theta/2 (cause spin two shape) with theta the argument of complex number ec
    return ec_mod, phi_e

#Plot the surface
def plot_halo(a, b, c) :
    x, y, z = halo(a, b, c)
    ax.set_box_aspect([a, b, c])
    ax.clear() #clearing previous plot
    ax.plot_surface(x, y, z,edgecolor='none', color='skyblue')
    plt.show()
#Plot the projection
def plot_proj(a, b, c, los) :
   ec_mod, phi_e = twoD_projection(a, b, c, los)
   plt.figure(figsize=(10, 6))
   plt.plot([abs(ec_mod)*.5*np.cos(phi_e), -abs(ec_mod)*.5*np.cos(phi_e)],
        [abs(ec_mod)*.5*np.sin(phi_e), -abs(ec_mod)*.5*np.sin(phi_e)], [0,0], linestyle = '-', marker = 'o', color = 'red')
   plt.title("2D projection of ellipsoid")
   plt.grid(True)
   plt.tight_layout()


#Window functions
def close(_) : #to close the window
  plt.close()

def update(_) : #updating display according to modifications
  # get values from cursors
  a=cursor_a.val
  b=cursor_b.val
  c=cursor_c.val
  plot_halo(a, b, c) #new plot

def reset(_) : #reset cursors
  cursor_a.reset()
  cursor_b.reset()
  cursor_c.reset()
  update(0) #updating display

def project(_) : #to have the 2D projection
    # get values from cursors
    los = los_button.value_selected
    a=cursor_a.val
    b=cursor_b.val
    c=cursor_c.val
    plot_proj(a, b, c, los) #plot of projection

def extract(_) :
   los = los_button.value_selected
   a=cursor_a.val
   b=cursor_b.val
   c=cursor_c.val
   eigenval, eigenvect = np.linalg.eigh(inertia_matrix(a, b, c))
   ec_mod, phi_e = twoD_projection(a, b, c, los)
   data = {
      'a' : [a],
      'b' : [b],
      'c' : [c],
      'los' : [los],
      'eigenvect' : [eigenvect.tolist()],
      'ec_mod' : [ec_mod],
      'phi_e' : [phi_e]
   }
   data_file = pd.DataFrame(data)
   data_file.to_csv(f'halo_{a:.2f}.{b:.2f}.{c:.2f}.{los}.csv', index=False) 

#Widget tracing
axe_a = plt.axes([0.1, 0.09, 0.65, 0.04]) #drawing crusor (with coordinates)
cursor_a= Slider(ax=axe_a, label='a', valmin=0, valmax=15, valinit=a0, track_color='darkgreen')
axe_b = plt.axes([0.1, 0.05, 0.65, 0.04]) #drawing crusor (with coordinates)
cursor_b= Slider(ax=axe_b, label='b', valmin=0, valmax=15, valinit=b0, track_color='darkgreen')
axe_c = plt.axes([0.1, 0.01, 0.65, 0.04]) #drawing crusor (with coordinates)
cursor_c= Slider(ax=axe_c, label='c', valmin=0, valmax=15, valinit=c0, track_color='darkgreen')
reset_canva=plt.axes([0.85, 0.05, 0.1, 0.03])
reset_button=Button(reset_canva,'Reset')
close_canva = plt.axes([0.85, 0.01, 0.1, 0.03])
close_button=Button(close_canva,'Close')
project_canva = plt.axes([0.84, 0.09, 0.12, 0.03])
project_button=Button(project_canva,'2D Projection')
los_canva = plt.axes([0.84, 0.15, 0.15, 0.15])
los_button=RadioButtons(los_canva,'xyz')
extract_canva = plt.axes([0.84, 0.55, 0.12, 0.03])
extract_button = Button(extract_canva, 'Extract halo')

#Widget activation
cursor_a.on_changed(update) # cursors
cursor_b.on_changed(update)
cursor_c.on_changed(update)
reset_button.on_clicked(reset) # reset
close_button.on_clicked(close) # close
project_button.on_clicked(project) # project
extract_button.on_clicked(extract) # extract halo data

#Initializing script
update(0)
plt.show(block=True)


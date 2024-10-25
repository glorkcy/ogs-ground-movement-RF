import numpy as np
from scipy.linalg import cholesky
import matplotlib.pyplot as plt
import scipy.spatial as sp
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import Mesh_settings 


# Import parameters from mesh_settings configuration file 

length_x = Mesh_settings.width_of_boundary # in m,   10
length_z = Mesh_settings.depth_of_boundary # in m,   60
npoint_x = Mesh_settings.horizontal_scale_of_fluctuation # number of cells,   10
npoint_z = Mesh_settings.vertical_scale_of_fluctuation # number of cells,   60

theta_x = Mesh_settings.horizontal_scale_of_fluctuation # 24.5
theta_z = Mesh_settings.vertical_scale_of_fluctuation  # 1.17

mean= Mesh_settings.stiffness_mean[0]  #1e7
slope = Mesh_settings.slope   # 2e5

trend_min = mean-(slope*(npoint_z-1)/2)

# Here we make the model
x = np.linspace(0., length_x, npoint_x).reshape(-1, 1)
z = np.linspace(0., length_z, npoint_z).reshape(-1, 1)

Dx = sp.distance_matrix(x, x)
Dz = sp.distance_matrix(z, z)

rho_corr_x = np.exp(-2 * Dx / theta_x)
rho_corr_z = np.exp(-2 * Dz / theta_z)

c_x = cholesky(rho_corr_x, lower=True)
c_z = cholesky(rho_corr_z, lower=True)

#U = np.random.lognormal(0, 1, size=[npoint_x, npoint_z])
U = np.random.normal(0, 1, size=[npoint_x, npoint_z])
#U = np.random.exponential(scale =1.0, size=[npoint_x, npoint_z])

X = np.dot(c_x, U.reshape(npoint_x,npoint_z)) 
X = X.reshape(npoint_x, npoint_z)
X = np.transpose(X) 

X = np.dot(c_z, X.reshape(npoint_z, npoint_x))
X = X.reshape(npoint_z, npoint_x)
X = np.transpose(X) 


def variance_depth_dependent (slope, z, trend_min):
    a = np.arange(0,z)
    y = 1+ a*slope/trend_min
    y = y[::-1]
    return y

slope_array = variance_depth_dependent (slope, npoint_z, trend_min)


# Multiplying arrays 
stiff_field = X * slope_array[np.newaxis, : ] *Mesh_settings.variance






#some notes (please skip)
#y = 1-x*slope/trend_max
#trend = mean + np.arange(-z/2, z/2)* slope
#y = 1+ slope /((mean-(slope*(z-1)/2)) +slope *x) 
#cov_min = Mesh_settings.variance /trend_min

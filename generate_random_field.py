import numpy as np
from scipy.linalg import cholesky
import scipy.spatial as sp
import Mesh_settings 


# Import parameters from mesh_settings configuration file 

length_x = Mesh_settings.width_of_boundary
length_z = Mesh_settings.depth_of_boundary


number_of_cells_x = npoint_x = Mesh_settings.number_of_cells_x
number_of_cells_z = npoint_z = Mesh_settings.number_of_cells_y

theta_x = Mesh_settings.horizontal_scale_of_fluctuation # 24.5
theta_z = Mesh_settings.vertical_scale_of_fluctuation  # 1.17

mean= Mesh_settings.stiffness_mean[0]  #1e7
slope = Mesh_settings.slope   # 2e5

number_of_layers = Mesh_settings.number_of_layers
COV = Mesh_settings.COV


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
np.random.seed(0) 
U = np.random.normal(0, 1, size=[npoint_x, npoint_z])

X = np.dot(c_x, U.reshape(npoint_x,npoint_z)) 
X = X.reshape(npoint_x, npoint_z)
X = np.transpose(X) 

X = np.dot(c_z, X.reshape(npoint_z, npoint_x))
X = X.reshape(npoint_z, npoint_x)
X = np.transpose(X) 


depth_of_layer =  length_z / number_of_layers

#slope, mean, hn, dl, vcpl
def depth_linear_dependent_layer (slope, m, dl, vcpl):
    x = np.arange(-vcpl/2, vcpl/2)
    delta = -slope * (dl / vcpl)
    y = (m + x * delta)
    return y

trend = depth_linear_dependent_layer (slope, mean, depth_of_layer, number_of_cells_z)

# Multiplying arrays 
stiff_field = X * trend[np.newaxis, : ] *COV


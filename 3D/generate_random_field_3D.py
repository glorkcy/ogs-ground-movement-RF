import numpy as np
from scipy.linalg import cholesky
import scipy.spatial as sp
import Mesh_settings_3D 


# Import parameters from Mesh_settings_3D_3D configuration file 

length_x = Mesh_settings_3D.length_of_boundary
length_y = Mesh_settings_3D.width_of_boundary
length_z = Mesh_settings_3D.depth_of_boundary


number_of_cells_x = npoint_x = Mesh_settings_3D.number_of_cells_x
number_of_cells_y = npoint_y = Mesh_settings_3D.number_of_cells_y
number_of_cells_z = npoint_z = Mesh_settings_3D.number_of_cells_z

theta_x = theta_y = Mesh_settings_3D.horizontal_scale_of_fluctuation # 24.5
theta_z = Mesh_settings_3D.vertical_scale_of_fluctuation  # 1.17

mean= Mesh_settings_3D.stiffness_mean[0]  #1e7
slope = Mesh_settings_3D.slope   # 2e5

number_of_layers = Mesh_settings_3D.number_of_layers
COV = Mesh_settings_3D.COV


# Here we make the model
x = np.linspace(0., length_x, npoint_x).reshape(-1, 1)
y = np.linspace(0., length_y, npoint_y).reshape(-1, 1)
z = np.linspace(0., length_z, npoint_z).reshape(-1, 1)

Dx = sp.distance_matrix(x, x)
Dy = sp.distance_matrix(y, y)
Dz = sp.distance_matrix(z, z)

rho_corr_x = np.exp(-2 * Dx / theta_x)
rho_corr_y = np.exp(-2 * Dy / theta_y)
rho_corr_z = np.exp(-2 * Dz / theta_z)

c_x = cholesky(rho_corr_x, lower=True)
c_y = cholesky(rho_corr_y, lower=True)
c_z = cholesky(rho_corr_z, lower=True)

#U = np.random.lognormal(0, 1, size=[npoint_x, npoint_z])
np.random.seed(0) 
U = np.random.normal(0, 1, size=[npoint_x, npoint_y, npoint_z])

# Matrix permutation according to Li et al. (2019)
X = np.dot(c_x, U.reshape(npoint_x, npoint_y * npoint_z))
X = X.reshape(npoint_x, npoint_y, npoint_z)
X = np.transpose(X, (1, 2, 0))

X = np.dot(c_y, X.reshape(npoint_y, npoint_z * npoint_x))
X = X.reshape(npoint_y, npoint_z, npoint_x)
X = np.transpose(X, (1, 2, 0))

X = np.dot(c_z, X.reshape(npoint_z, npoint_x * npoint_y))
X = X.reshape(npoint_z, npoint_x, npoint_y)
#X = np.transpose(X, (1, 2, 0))

depth_of_layer =  length_z / number_of_layers

#slope, mean, hn, dl, vcpl
def depth_linear_dependent_layer (slope, m, dl, vcpl):
    x = np.arange(-vcpl/2, vcpl/2)
    delta = -slope * (dl / vcpl)
    y = (m + x * delta)
    return y

trend = depth_linear_dependent_layer (slope, mean, depth_of_layer, number_of_cells_z)

# make the trend to a 3Darray
trend = trend[:,np.newaxis]
trend = trend[:,np.newaxis]

# Multiplying arrays 
stiff_field = X * trend  *COV


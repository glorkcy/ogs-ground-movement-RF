import numpy as np
from scipy.linalg import cholesky
import scipy.spatial as sp
import Mesh_settings


# Import parameters from Mesh_settings

length_x = length_y = Mesh_settings.width_of_boundary
length_z = Mesh_settings.depth_of_boundary

number_of_cells_x = number_of_cells_y = npoint_x = npoint_y =Mesh_settings.number_of_cells_x
number_of_cells_z = npoint_z = Mesh_settings.number_of_cells_y

theta_x = theta_y = Mesh_settings.horizontal_scale_of_fluctuation # 24.5
theta_z = Mesh_settings.vertical_scale_of_fluctuation  # 1.17

mean= Mesh_settings.stiffness_mean[0]  #1e7
slope = Mesh_settings.slope   # 2e5

number_of_layers = Mesh_settings.number_of_layers
COV = Mesh_settings.COV

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
np.random.seed(1) 
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
X = np.transpose(X, (1, 2, 0))

X = X[:,0,:]

depth_of_layer =  length_z / number_of_layers

#slope, mean, hn, dl, vcpl
def depth_linear_dependent_layer (slope, m, dl, vcpl):
    x = np.arange(-vcpl/2, vcpl/2)
    delta = -slope * (dl / vcpl)
    y = (m + x * delta)
    return y

def formula_derived_layer (depth_of_layer, number_of_cells_z): #Grundbau Taschenbuch Tabelle 1
    lower_limit_sigma = 0.3/(1-0.3)*(1-0.35)*2600*9.81*1 #0.3/(1-0.3)*
    upper_limit_sigma = 0.3/(1-0.3)*(1-0.35)*2600*9.81*depth_of_layer
    
    v = 225
    w = 0.675
    
    total_stress = np.linspace(upper_limit_sigma, lower_limit_sigma, number_of_cells_z, endpoint=True, retstep=False, dtype=None)

    E = v * 100000 * (total_stress/100000)**w

    return E

trend = formula_derived_layer(depth_of_layer, number_of_cells_z)

# Multiplying arrays 
stiff_field = X * trend[np.newaxis, : ] *COV


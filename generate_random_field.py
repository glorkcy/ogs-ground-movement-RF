import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as sp
from _out import Mesh_settings

# Parameters
lx, ly, lz = 10, 12, 60  # Length of soil in x, y, z directions (meters)   #100,100,60
ex, ey, ez = 10, 12, 60  # Number of elements in each direction             #100,100,600
nx, ny, nz = ex + 1, ey + 1, ez + 1  # Number of nodes in each direction
scales = (22, 22, 1)  # Correlation scales in x, y, z directions (meters)
cov = 0  # Coefficient of variation  #0.35
layers =  []
initial_GWT = - Mesh_settings.initial_GWT_depth

# # Function to compute depth-dependent stiffness modulus (linear relationship)
# stiffness_mean = 1e7  # Mean stiffness modulus (Pa)
# slope = 2e5  # Stiffness gradient per meter
# def depth_linear_dependent_layer(stiffness_mean, slope, lz, ez):
#     z_positions = np.arange(-ez / 2, ez / 2)  # Depth positions centered around zero
#     delta = slope * (lz / ez)
#     E = stiffness_mean + z_positions * delta
#     return E

def formula_derived_layer(lz, nz, layers, initial_GWT):
    """
    Compute the stiffness modulus (E) with multiple layers, each having different v and w.

    Parameters:
    - lz (float): Total depth
    - nz (int): Number of depth points
    - layers (list of tuples): Each tuple contains (depth_start, depth_end, v, w, gamma)
    - unit_weight (float): Default unit weight of soil (2600 kg/m³ * 9.81 m/s²)

    Returns:
    - E (numpy array): Stiffness modulus for each depth
    """
    depth_intervals = np.linspace(0, lz, nz)  # Generate depth positions
    total_stress = np.zeros(nz)  # Initialize stress array
    effective_stress = np.zeros(nz)  # Initialize stress array
    E = np.zeros(nz)  # Initialize stiffness modulus array
    
    # Compute total stress at each depth considering soil layers
    for i, depth in enumerate(depth_intervals):
        stress = 0  # Initialize pore pressure

        for (depth_start, depth_end, v, w, density, nu, phi) in layers:
            if depth >= depth_start:
                # Contribution of this layer to stress
                dz = min(depth, depth_end) - depth_start  # Depth within this layer
                stress += dz * density * 9.81 * (1 - phi)  # Accumulate total weight

            # Check if current depth is below water table 
        # Calculate pore pressure if below water table
        
        total_stress[i] = stress
        initial_GWT_depth = -initial_GWT
        if depth >= initial_GWT_depth:
            water_height = depth - initial_GWT_depth
            pore_pressure = water_height * 9.81 * 1000  # Pa
            effective_stress[i] = total_stress[i] - pore_pressure
        else:
            effective_stress[i] = total_stress[i]  # Above water table

    # Compute stiffness modulus E using v and w for each layer
    for (depth_start, depth_end, v, w, density, nu, phi) in layers:
        mask = (depth_intervals >= depth_start) & (depth_intervals <= depth_end)  #mask is a Boolean array where True means that the corresponding depth_intervals fall within the current layer
        E[mask] = (v * 1e5 * (effective_stress[mask] / 1e5) ** w ) * (1+nu)*(1-2*nu)/(1-nu)  # Compute E for this layer and young modulus conversion

    return E


# Stiffness modulus trend along the z-direction
#stiffness_modulus_z = depth_linear_dependent_layer(stiffness_mean, slope, lz, ez)
stiffness_modulus_z = formula_derived_layer(lz, nz, layers,initial_GWT)

def create_1d_correlation_matrix(nodes, scale, length, acf_type='-'):
    """
    Generates 1D correlation matrix for different ACF models.

    Args:
        nodes: Number of nodes.
        scale: Scale of fluctuation
        length: Length of the domain.
        acf_type: Type of ACF ('exponential', 'gaussian', 'markov', 'binary_noise'). 
                 Defaults to 'exponential'.

    Returns:
        A numpy array representing the correlation matrix.
    """
    indices = np.linspace(0., length, nodes).reshape(-1, 1)
    distances = sp.distance_matrix(indices, indices)

    if acf_type == 'exponential':
        ACF = np.exp(-2 * np.abs(distances) / scale)
    elif acf_type == 'gaussian':
        ACF = np.exp(- np.pi * distances**2 / scale**2)
    elif acf_type == 'second_markov':
        ACF = (1 + 4 * np.abs(distances) / scale) * np.exp(-4 * np.abs(distances) / scale)
    elif acf_type == 'binary_noise':
        ACF = np.where(distances <= scale, 1 - distances / scale, 0)
    elif acf_type == 'cosine_exponential':
        ACF = np.exp(-np.abs(distances)/ scale) * np.cos (np.abs(distances)/ scale)
    elif acf_type == 'spherical':
        ACF = np.where(distances <= 4/3*scale, 1 -9/8* np.abs(distances / scale) +27/128 *(np.abs(distances/scale)**3), 0)
    else:
        raise ValueError(f"Invalid ACF type: {acf_type}")

    if acf_type == 'gaussian':  # Ensure positive definiteness for Gaussian
        eigenvalues = np.linalg.eigvals(ACF)
        if np.any(eigenvalues <= 0):
            ACF += np.eye(ACF.shape[0]) * 1e-14

    return ACF

def stepwise_cmd_1d_random_field(nz, lz, scale_z, cov, stiffness_modulus_z, acf_type='-', distribution='-'):
    """
    Generate a 1D random field using the stepwise CMD method.
    """

    # Create 1D correlation matrix
    Rz = create_1d_correlation_matrix(nz, scale_z, lz, acf_type)
    
    # Perform Cholesky decomposition
    Lz = np.linalg.cholesky(Rz)
    
    # Generate standard normal random values
    np.random.seed(12517)
    U = np.random.normal(size=(nz,))
    
    # Apply Cholesky decomposition
    X = Lz @ U

    # Apply mean and standard deviation along depth
    if distribution == 'lognormal':
        mean_z = np.log(stiffness_modulus_z) - 0.5 * np.log(1 + cov**2)
        std_dev_z = np.sqrt(np.log(1 + cov**2)) * np.ones_like(mean_z)
    elif distribution == 'normal':
        mean_z = stiffness_modulus_z
        std_dev_z = cov * stiffness_modulus_z
    else:
        raise ValueError("Invalid distribution. Choose 'lognormal' or 'normal'.")
    
    X = mean_z + std_dev_z * X
    
    # Convert to lognormal if needed
    if distribution == 'lognormal':
        X = np.exp(X)
    
    return X


# Stepwise CMD method for random field generation
def stepwise_cmd_3d_random_field(nx, ny, nz, lx, ly, lz, scales,cov, stiffness_modulus_z, acf_type='-', distribution='-'):
    """
    Generate a 3D random field using the stepwise CMD method.
    """
    scale_x, scale_y, scale_z = scales

    Rx = create_1d_correlation_matrix(nx, scale_x, lx, acf_type)
    Ry = create_1d_correlation_matrix(ny, scale_y, ly, acf_type)
    Rz = create_1d_correlation_matrix(nz, scale_z, lz, acf_type)
    
    # Perform Cholesky decomposition
    Lx = np.linalg.cholesky(Rx) #+ 1e-6 * np.eye(nx)
    Ly = np.linalg.cholesky(Ry)
    Lz = np.linalg.cholesky(Rz)

    # Generate standard normal random values
    np.random.seed(123)
    U = np.random.normal(size=(nx, ny, nz))

    # Apply Cholesky decomposition stepwise
    X = np.reshape(U, (nx, -1))
    X = Lx @ X
    X = np.reshape(X, (nx, ny, nz))
    X = np.transpose(X, (1, 2, 0))

    X = np.reshape(X, (ny, -1))
    X = Ly @ X
    X = np.reshape(X, (ny, nz, nx))
    X = np.transpose(X, (1, 2, 0))

    X = np.reshape(X, (nz, -1))
    X = Lz @ X
    X = np.reshape(X, (nz, nx, ny))
    X = np.transpose(X, (1, 2, 0))

    # Apply mean and standard deviation along depth
    if distribution == 'lognormal':
        mean_z = np.log(stiffness_modulus_z) - 0.5 * np.log(1 + cov**2)
        std_dev_z = np.sqrt(np.log(1 + cov**2)) * np.ones_like(mean_z)
    elif distribution == 'normal':
        mean_z = stiffness_modulus_z
        std_dev_z = cov * stiffness_modulus_z
    else:
        raise ValueError("Invalid distribution. Choose 'lognormal' or 'normal'.")

    for k in range(nz):
        X[:, :, k] = mean_z[k] + std_dev_z[k] * X[:, :, k]

    # Convert to lognormal if needed
    if distribution == 'lognormal':
        X = np.exp(X)

    return X

# Generate fields
#node_field_gaussian_lognormal = stepwise_cmd_3d_random_field(nx, ny, nz, lx, ly, lz, scales, cov, stiffness_modulus_z, acf_type = 'gaussian', distribution = 'lognormal')
#node_field_exponential_lognormal = stepwise_cmd_3d_random_field(nx, ny, nz, lx, ly, lz, scales,cov, stiffness_modulus_z, acf_type='Exponential', distribution='lognormal')
#node_field_gaussian_normal = stepwise_cmd_3d_random_field(nx, ny, nz, lx, ly, lz, scales,cov, stiffness_modulus_z, acf_type='Gaussian', distribution='normal')
#node_field_exponential_normal = stepwise_cmd_3d_random_field(nx, ny, nz, lx, ly, lz, scales, cov, stiffness_modulus_z, acf_type='Exponential', distribution='normal')

def convert_to_element_field_1d(node_field, ez):
    """
    Convert a node-based random field to an element-based field by averaging.
    """
    element_field = np.zeros(ez)

    for i in range(ez):
        # Average values of the 8 corner nodes for the element
        element_field[i] = np.mean(node_field[i:i+2])
    
    return element_field


# Convert node-based fields to element-based fields
def convert_to_element_field(node_field, ex, ey, ez):
    """
    Convert a node-based random field to an element-based field by averaging.
    """
    element_field = np.zeros((ex, ey, ez))

    for i in range(ex):
        for j in range(ey):
            for k in range(ez):
                # Average values of the 8 corner nodes for the element
                element_field[i, j, k] = np.mean(node_field[i:i+2, j:j+2, k:k+2])

    return element_field


#element_field_gaussian_lognormal = convert_to_element_field(node_field_gaussian_lognormal, ex, ey, ez)
#element_field_exponential_lognormal = convert_to_element_field(node_field_exponential_lognormal, ex, ey, ez)
#element_field_gaussian_normal = convert_to_element_field(node_field_gaussian_normal, ex, ey, ez)
#element_field_exponential_normal = convert_to_element_field(node_field_exponential_normal, ex, ey, ez)


#selected_field = element_field_exponential_lognormal[:,0,:]
#stiff_field=np.flip(selected_field, axis=1)


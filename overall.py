import os
from pathlib import Path
import pyvista as pv
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk
import ogs as ogs
import xml.etree.ElementTree as ET
import ogstools as ot
import vtuIO
import matplotlib.pyplot as plt
import scipy.spatial as sp
import scipy.optimize as opt
from scipy.optimize import curve_fit, minimize, differential_evolution

class MeshGeneratorWithStiffness:
    def __init__(self, out_dir="_out", mesh_settings=None):
        """
        Initialize the mesh generator with stiffness properties.
        
        Args:
            out_dir (str): Output directory for the generated mesh files
            mesh_settings (module): Module containing mesh settings (lx, ly, lz, ex, ey, ez, layers, initial_GWT_depth, duration)
                                 If None, will try to import Mesh_settings module
        """
        self.out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", out_dir))
        self.mesh_settings = mesh_settings
        self.mesh = None
        self.reader_data = None
        
        # Default parameters
        self.scales = (22, 1.17)  # Correlational length (horizontal, vertical)
        self.acf_type = "exponential"
        self.distribution = "lognormal"
        self.cov = 0.3
        self.g = -9.81  # gravity
        self.output_file = "domain_with_stiffness.vtu"
        
    def setup_environment(self):
        """Set up the output directory and working environment."""

        os.chdir(self.out_dir)
        import Mesh_settings
        self.mesh_settings = Mesh_settings

    def generate_empty_mesh(self):
        """Generate the initial empty mesh structure."""
        lx, ly, lz = self.mesh_settings.lx, self.mesh_settings.ly, self.mesh_settings.lz
        ex, ey, ez = self.mesh_settings.ex, self.mesh_settings.ey, self.mesh_settings.ez
        
        file = 'empty_mesh.vtu'
        
        if ey > 1:  # 3D case
            ogs.cli.generateStructuredMesh(o=file, e='hex', lx=lx, ly=ly, lz=lz, nx=ex, ny=ey, nz=ez, oz=-lz)
        else:  # 2D case
            ogs.cli.generateStructuredMesh(o=file, e='quad', lx=lx, ly=lz, nx=ex, ny=ez, oy=-lz)
        
        # Identify subdomains
        ogs.cli.identifySubdomains(m=file, s='1e-6', f=file)
        
        # PyVista settings
        pv.set_plot_theme("document")
        pv.set_jupyter_backend("static")
        
        # Read the empty mesh
        self.mesh = pv.read(file)
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(file)
        reader.Update()
        self.reader_data = reader.GetOutput()
        
    def formula_derived_layer(self, lz, nz, layers, initial_GWT):
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

    def create_1d_correlation_matrix(self, nodes, scale, length, acf_type='-'):
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

    def stepwise_cmd_1d_random_field(self, nz, lz, scale_z, cov, stiffness_modulus_z, acf_type='-', distribution='-'):
        """
        Generate a 1D random field using the stepwise CMD method.
        """

        # Create 1D correlation matrix
        Rz = self.create_1d_correlation_matrix(nz, scale_z, lz, acf_type)
        
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

    def stepwise_cmd_3d_random_field(self, nx, ny, nz, lx, ly, lz, scales,cov, stiffness_modulus_z, acf_type='-', distribution='-'):
        """
        Generate a 3D random field using the stepwise CMD method.
        """
        scale_x, scale_y, scale_z = scales

        Rx = self.create_1d_correlation_matrix(nx, scale_x, lx, acf_type)
        Ry = self.create_1d_correlation_matrix(ny, scale_y, ly, acf_type)
        Rz = self.create_1d_correlation_matrix(nz, scale_z, lz, acf_type)
        
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

    def convert_to_element_field_1d(self, node_field, ez):
        """
        Convert a node-based random field to an element-based field by averaging.
        """
        element_field = np.zeros(ez)

        for i in range(ez):
            # Average values of the 8 corner nodes for the element
            element_field[i] = np.mean(node_field[i:i+2])
        
        return element_field

    def convert_to_element_field(self, node_field, ex, ey, ez):
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

    def generate_stiffness_field(self):
        """Generate the random stiffness field."""
        lx, ly, lz = self.mesh_settings.lx, self.mesh_settings.ly, self.mesh_settings.lz
        ex, ey, ez = self.mesh_settings.ex, self.mesh_settings.ey, self.mesh_settings.ez
        nx, ny, nz = ex + 1, ey + 1, ez + 1
        
        trend = self.formula_derived_layer(
            lz, nz, self.mesh_settings.layers, -self.mesh_settings.initial_GWT_depth
        )
        
        if ex == 1 and ey == 1:  # 1D case
            field = self.stepwise_cmd_1d_random_field(
                nz, lz,
                self.scales[1],  # Correlation scales
                self.cov,
                trend,
                acf_type=self.acf_type,
                distribution=self.distribution
            )
            stiff_field = self.convert_to_element_field_1d(field, ez)
            stiff_field = np.flip(stiff_field, axis=0)
        else:
            field = self.stepwise_cmd_3d_random_field(
                nx, ny, nz,
                lx, ly, lz,
                (self.scales[0], self.scales[0], self.scales[1]),
                self.cov,
                trend,
                acf_type=self.acf_type,
                distribution=self.distribution
            )
            
            element_field = self.convert_to_element_field(field, ex, ey, ez)
            
            if ey == 1:
                stiff_field = element_field[:, 0, :]
                stiff_field = np.flip(stiff_field, axis=1)
            else:
                stiff_field = element_field[:, :, :]
                stiff_field = np.flip(stiff_field, axis=2)
        
        return stiff_field
    
    def calculate_material_properties(self):
        """Calculate material properties based on depth and layers."""
        lx, ly, lz = self.mesh_settings.lx, self.mesh_settings.ly, self.mesh_settings.lz
        ex, ey, ez = self.mesh_settings.ex, self.mesh_settings.ey, self.mesh_settings.ez
        initial_GWT = -self.mesh_settings.initial_GWT_depth
        layers = self.mesh_settings.layers
        
        dz = lz / ez
        depths = np.linspace(-lz + dz/2, -dz/2, ez)  # cell centers
        
        # Initialize arrays
        density_1D_array = np.zeros(ez)
        nu_1D_array = np.zeros(ez)
        phi_1D_array = np.zeros(ez)
        initial_pore_pressure = np.zeros(ez)
        initial_stress_xx_1D_array = np.zeros(ez)
        initial_stress_yy_1D_array = np.zeros(ez)
        
        for i, depth in enumerate(depths):
            sigma_yy = 0
            nu_at_depth = 0
            density_at_depth = 0
            phi_at_depth = 0

            for z_start, z_end, v, w, density, nu, phi in layers:
                layer_top = -z_start
                layer_bottom = -z_end

                if layer_bottom >= depth:
                    thickness = layer_top - layer_bottom
                    if layer_top >= initial_GWT and layer_bottom >= initial_GWT:
                        gamma = (1 - phi) * density * self.g
                    elif layer_top < initial_GWT and layer_bottom < initial_GWT:
                        gamma = ((1 - phi) * density + phi * 1000) * self.g
                    else:
                        above_thickness = layer_top - initial_GWT
                        below_thickness = initial_GWT - layer_bottom
                        gamma_above = (1 - phi) * density * self.g
                        gamma_below = ((1 - phi) * density + phi * 1000) * self.g
                        sigma_yy += gamma_above * above_thickness
                        sigma_yy += gamma_below * below_thickness
                        continue
                    sigma_yy += gamma * thickness

                elif layer_bottom < depth <= layer_top:
                    thickness = layer_top - depth
                    nu_at_depth = nu
                    density_at_depth = density
                    phi_at_depth = phi
                    if layer_top >= initial_GWT and depth >= initial_GWT:
                        gamma = (1 - phi) * density * self.g
                    elif layer_top < initial_GWT and depth < initial_GWT:
                        gamma = ((1 - phi) * density + phi * 1000) * self.g
                    else:
                        above_thickness = layer_top - initial_GWT
                        below_thickness = initial_GWT - depth
                        gamma_above = (1 - phi) * density * self.g
                        gamma_below = ((1 - phi) * density + phi * 1000) * self.g
                        sigma_yy += gamma_above * above_thickness
                        sigma_yy += gamma_below * below_thickness
                        thickness = 0
                    if thickness > 0:
                        sigma_yy += gamma * thickness

            # Store material properties
            density_1D_array[i] = density_at_depth
            nu_1D_array[i] = nu_at_depth
            phi_1D_array[i] = phi_at_depth

            # Compute stress components
            sigma_xx = nu_at_depth / (1 - nu_at_depth) * sigma_yy
            initial_stress_yy_1D_array[i] = sigma_yy
            initial_stress_xx_1D_array[i] = sigma_xx

            # Compute pore pressure
            if depth < initial_GWT:
                initial_pore_pressure[i] = 1000 * (initial_GWT - depth) * self.g
            else:
                initial_pore_pressure[i] = 0.0
        
        return (
            density_1D_array, nu_1D_array, phi_1D_array,
            initial_pore_pressure, initial_stress_xx_1D_array, initial_stress_yy_1D_array
        )
    
    def add_parameters_to_mesh(self):
        """Add all calculated parameters to the mesh."""
        ex, ey, ez = self.mesh_settings.ex, self.mesh_settings.ey, self.mesh_settings.ez
        
        # Generate and add stiffness field
        stiff_field = self.generate_stiffness_field()
        stiffness_1D_array = stiff_field.flatten('F')
        stiffness_vtk = numpy_to_vtk(stiffness_1D_array)
        stiffness_vtk.SetName('Stiffness')
        self.reader_data.GetCellData().AddArray(stiffness_vtk)
        
        # Calculate material properties
        (
            density_1D_array, nu_1D_array, phi_1D_array,
            initial_pore_pressure, initial_stress_xx_1D_array, initial_stress_yy_1D_array
        ) = self.calculate_material_properties()
        
        # Flatten and repeat arrays to match mesh shape
        n_cells = ex * ey * ez
        density_1D_array = np.repeat(density_1D_array, ex * ey)
        nu_1D_array = np.repeat(nu_1D_array, ex * ey)
        phi_1D_array = np.repeat(phi_1D_array, ex * ey)
        initial_pore_pressure = np.repeat(initial_pore_pressure, ex * ey)
        initial_stress_xx_1D_array = np.repeat(initial_stress_xx_1D_array, ex * ey)
        initial_stress_yy_1D_array = np.repeat(initial_stress_yy_1D_array, ex * ey)
        
        # Combine stress components
        if ey == 1:
            initial_stress_1D_array = np.column_stack((
                initial_stress_xx_1D_array,
                initial_stress_yy_1D_array,
                initial_stress_xx_1D_array,
                np.zeros_like(initial_stress_xx_1D_array)
            ))
        else:
            initial_stress_1D_array = np.column_stack((
                initial_stress_xx_1D_array,
                initial_stress_xx_1D_array,
                initial_stress_yy_1D_array,
                np.zeros_like(initial_stress_xx_1D_array),
                np.zeros_like(initial_stress_xx_1D_array),
                np.zeros_like(initial_stress_xx_1D_array)
            ))
        
        # Add arrays to VTK data
        density_vtk = numpy_to_vtk(density_1D_array)
        nu_vtk = numpy_to_vtk(nu_1D_array)
        phi_vtk = numpy_to_vtk(phi_1D_array)
        initial_stress_vtk = numpy_to_vtk(initial_stress_1D_array)
        
        density_vtk.SetName('Density')
        nu_vtk.SetName('nu')
        phi_vtk.SetName('Porosity')
        initial_stress_vtk.SetName('initial_stress')
        
        self.reader_data.GetCellData().AddArray(density_vtk)
        self.reader_data.GetCellData().AddArray(nu_vtk)
        self.reader_data.GetCellData().AddArray(phi_vtk)
        self.reader_data.GetCellData().AddArray(initial_stress_vtk)
    
    def save_mesh(self):
        """Save the final mesh with all parameters."""
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputData(self.reader_data)
        writer.SetFileName(self.output_file)
        writer.Write()
        return str(self.out_dir / self.output_file)
    
    def generate(self):
        """Main method to generate the complete mesh with all properties."""
        self.setup_environment()
        self.generate_empty_mesh()
        self.add_parameters_to_mesh()
        return self.save_mesh()

class OGSModelRunner:
    def __init__(self, out_dir=None):
        """
        Initialize the OGS model runner.
        
        Args:
            out_dir (str, optional): Output directory path. Defaults to "_out" or environment variable.
        """
        self.out_dir = Path(out_dir if out_dir else os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
        self.ey = None
        self.duration = None
        self.delta_time = None
        self.GWT = None
        self.lx = None
        self.ly = None
        self.lz = None
        
    def load_settings(self):
        """Load settings from Mesh_settings module."""
        import Mesh_settings
        self.ey = Mesh_settings.ey
        self.duration = Mesh_settings.duration
        self.delta_time = Mesh_settings.delta_time
        self.GWT = -Mesh_settings.initial_GWT_depth
        self.lx = Mesh_settings.lx
        self.ly = Mesh_settings.ly
        self.lz = Mesh_settings.lz
        
    def modify_geometry(self):
        """
        Modify the geometry file (2D or 3D) based on the settings.
        """
        self.load_settings()
        
        if self.ey == 1:  # 2D case
            self._modify_2d_geometry()
        else:  # 3D case
            self._modify_3d_geometry()
    
    def _modify_2d_geometry(self):
        """Modify the 2D rectangle geometry file."""
        tree = ET.parse('rectangle.gml')
        root = tree.getroot()
        
        points = root.findall('./points/point')
        for p in points:
            if p.get('id') in {'1', '2'}:
                p.set('x', str(self.lx))
            if p.get('id') in {'2', '3'}:
                p.set('y', str(-self.lz))
                
        tree.write('rectangle.gml')
    
    def _modify_3d_geometry(self):
        """Modify the 3D cuboid geometry file."""
        tree = ET.parse('cuboid.gml')
        root = tree.getroot()
        
        points = root.findall('./points/point')
        for p in points:
            if p.get('id') in {'4', '5', '6', '7'}:
                p.set('x', str(self.lx))
            if p.get('id') in {'2', '3', '6', '7'}:
                p.set('y', str(self.ly))    
            if p.get('id') in {'0', '3', '4', '7'}:
                p.set('z', str(-self.lz))
                
        tree.write('cuboid.gml')
    
    def configure_and_run_model(self):
        """
        Configure the OGS model with the current settings and run it.
        """
        # Determine the project file based on dimensionality
        project_file = "input_2D.prj" if self.ey == 1 else "input_3D.prj"
        
        # Initialize OGS model
        model = ot.Project(input_file=project_file, output_file=project_file)
        
        # Modify groundwater table expression
        gwt_expression = f'if (y >= {self.GWT}, 0, 1000*(-y {self.GWT})*9.81)'
        model.replace_text(
            gwt_expression,
            xpath="./parameters/parameter/expression",
            occurrence=0
        )
        
        # Modify timesteps
        model.replace_text(
            str(int(self.duration / self.delta_time)),
            xpath='./time_loop/processes/process/time_stepping/t_end'
        )
        model.replace_text(
            str(int(self.duration / self.delta_time)),
            xpath='./time_loop/processes/process/time_stepping/timesteps/pair/repeat'
        )
        
        # Write input and run model
        model.write_input()
        model.run_model(logfile="out.txt")
    
    def run(self):
        """
        Execute the complete workflow:
        1. Change to output directory
        2. Modify geometry
        3. Configure and run model
        4. Return to original directory
        """
        original_dir = os.getcwd()
        try:
            self.modify_geometry()
            self.configure_and_run_model()
        finally:
            os.chdir(original_dir)

class GroundwaterUpliftAnalyzer:
    def __init__(self):
        self.ey = None
        self.ez = None
        self.lz = None
        self.initial_GWT_depth = None
        self.duration = None
        self.delta_time = None
        self.layers = None
        self.dim = None
        self.ax = None
        self.points = None
        self.sigma_ref = 100000.0  # Example reference stress in Pascals
        self.g = 9.81  # Gravity (m/s²)
        self.pvdfile = None
        self.accumulated_u = None
        self.total_u_at_timestep = None
        self.particular_u_at_timestep_half_saturated = None
        self.particular_u_at_timestep_fully_saturated = None
        self.stiffness_modulus = None
        
    def load_settings(self):
        import Mesh_settings
        
        # Get some data
        self.ey = Mesh_settings.ey
        self.ez = Mesh_settings.ez
        self.lz = Mesh_settings.lz
        self.initial_GWT_depth = Mesh_settings.initial_GWT_depth
        self.duration = Mesh_settings.duration
        self.delta_time = Mesh_settings.delta_time
        self.layers = Mesh_settings.layers

    def validate_model(self):
        # Now we just keep it as one layer model first. Later the script should be edited to fit multi-layers
        if len(self.layers) != 1:
            raise ValueError('It only applies to 1 layer model.')
            
        self.z_start, self.z_end, self.v_original, self.w_original, self.rho, self.nu, self.porosity = self.layers[0]

    def setup_coordinates(self):
        if self.ey == 1:  # 2D
            self.dim = 2                       # Dimension
            self.ax = 1                        # Axis
            self.points = {"pt0": (0, 0)}      # Define point of interest (x,z), z has to be 0
        else:  # 3D
            self.dim = 3                    
            self.ax = 2                           
            self.points = {"pt0": (0, 0, 0)}  # (x,y,z)

    def compute_depths(self):
        # Compute the depth affected by the rising of groundwater
        depth_upper = self.duration - self.initial_GWT_depth  
        depth_lower = -self.initial_GWT_depth 
        resolution_mesh = self.lz/self.ez
        self.depth = np.arange(depth_upper, depth_lower, -self.delta_time)
        self.depth_original = np.arange(depth_upper, depth_lower, -resolution_mesh)

    def load_displacement_data(self):
        # Read the pvd file
        self.pvdfile = vtuIO.PVDIO("uplift.pvd", dim=self.dim, interpolation_backend="scipy")
        np.set_printoptions(precision=12)  # Set numerical precision of data
        
        # Read the accumulated displacement data
        self.accumulated_u = self.pvdfile.read_time_series("displacement", self.points, 
                                                         interpolation_method="nearest")

    def compute_displacements(self):
        # Calculate the time-specific displacement of all cells
        self.total_u_at_timestep = (self.accumulated_u["pt0"][1:, self.ax] - 
                                   self.accumulated_u["pt0"][:-1, self.ax])

        # Calculate the displacement of each GW-uplift-experiencing cell
        self.particular_u_at_timestep_half_saturated = self._compute_pts(self.total_u_at_timestep)
        self.particular_u_at_timestep_fully_saturated = 2 * self.particular_u_at_timestep_half_saturated

    def _compute_pts(self, total_u_at_timestep):
        n = len(total_u_at_timestep)
        pts = np.zeros(n)
        for i in range(n):
            pts[i] = total_u_at_timestep[i] - 2 * np.sum(pts[:i])  # minus the displacement of all bottom fully-saturated cells
        return pts

    def calculate_stiffness_modulus(self):
        # Calculate Es: mesh height^2 * unit weight of water * (1-porosity) / (2* displacement)
        self.stiffness_modulus = (self.delta_time**2 * 9.81e3 * (1-self.porosity) / 
                                (self.particular_u_at_timestep_fully_saturated))
        self.stiffness_modulus = self.stiffness_modulus[::-1]

    def run_analysis(self):
        self.load_settings()
        self.validate_model()
        self.setup_coordinates()
        self.compute_depths()
        self.load_displacement_data()
        self.compute_displacements()
        self.calculate_stiffness_modulus()
        #return self.stiffness_modulus
        return {
            'stiffness_modulus': self.stiffness_modulus,
            'particular_u_at_timestep_half_saturated': self.particular_u_at_timestep_half_saturated
        }

def main():
    # Step 1: Generate mesh with stiffness
    mesh_generator = MeshGeneratorWithStiffness()
    mesh_file = mesh_generator.generate()
    print(f"Generated mesh file: {mesh_file}")
    
    # Step 2: Run OGS model
    ogs_runner = OGSModelRunner()
    ogs_runner.run()
    print("OGS model run completed")
    
    # Step 3: Analyze results and get stiffness
    analyzer = GroundwaterUpliftAnalyzer()
    results = analyzer.run_analysis()
    print("Calculated stiffness modulus:", results['stiffness_modulus'])
    print("Individual displacement of GW-experiencing half-saturated cell :", results['particular_u_at_timestep_half_saturated'])

if __name__ == "__main__":
    main()
    

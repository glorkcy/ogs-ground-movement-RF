Simulation of the uplift of ground, due to groundwater lifting.

### Brief Introduction:
In this example, groundwater table rises from -45m to -5m (0.5m per timestep).
The output graph is uplift displacement, which can be used to back estimate to the original stiffness equation.

### Available Features
- 2D or 3D
- Soil layers
- Gaussian / Exponential autocorrelation function
- Lognormal / Normal probability density function
- Initial groundwater table

### For Back-estimation
- Please keep the setting to 1D, 1-layer, Groundwater table at the bottommost

### Steps
- Step 1: Edit Mesh_settings.py inside the _out folder.
- Step 2: Run generate_mesh_cell. 
- Step 3: Run run_ogs.
- Step 4: Run run_pvd, to see the displacement graphs.

### Prerequisite: 
pip install ogs ogstools VTUinterface pyvista vtk 

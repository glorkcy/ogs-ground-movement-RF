This is a simple simulation of the effect of uprising groundwater table to the uplift of ground. 

### Brief Introduction:
The 2D mesh has size of 10m (horizontal) and 60m (vertically underground). 
Groundwater table is assumed to rise 1m per timestep.
There are in total 45 timesteps, so the water rises from -50m to -5m. 

### To understand the effect of random field: 

Step 1: Edit generate_random_field.py to input a different random field model.
(a numpy array, with shape: number of horizontal cells x number of vertical cells)

Step 2: Run overall.ipynb to generate uplift graphs

Prerequisite: pip install ogs6py VTUinterface nbformat nbconvert pyvista vtk 

P.S. All the parameters (e.g. stiffness mean, variance, ground properties etc.) are subject to change. In the meantime, a workflow is designed first.


import gstools as gs
import Mesh_settings 

# Import parameters from mesh_settings configuration file 

x = range(Mesh_settings.number_of_cells_x)
y = range(Mesh_settings.number_of_cells_y)
variance = Mesh_settings.variance     # for making random field
h_cl = Mesh_settings.horizontal_scale_of_fluctuation 
v_cl = Mesh_settings.vertical_scale_of_fluctuation 
random_seed = 12345 # just select a random number

# Here we make the model

stiff_model = gs.Exponential(dim=2, var= variance, len_scale=[h_cl,v_cl])
srf = gs.SRF(stiff_model, seed=random_seed)
stiff_field = srf.structured([x,y])

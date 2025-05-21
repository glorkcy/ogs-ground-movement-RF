""" lx, ly, lz are the mesh dimension (in meter), while nx, ny, nz are the number of elements.
    Setting ex and ey = 1, model is 1D, setting just ey = 1, model is 2D.
    Vertical correlational length is roughly 1m now, so we would set the vertical mesh size as 0.5m, to keep it smaller than the correlational length."""

lx = 1 
ly = 1
lz = 45 
ex = 1
ey = 1
ez = 90

# this annotate z_start, z_end, v, w, density, nu, porosity
layers = [
    (0, 45, 225, 0.675, 1750, 0.35, 0.4) 
]

initial_GWT_depth = 45
duration = 40
delta_time = 0.5

# other possible layer options
#(20, 40, 400, 0.6, 2000, 0.27, 0.33),  
#(40, 61, 600, 0.5, 2200, 0.2, 0.25) 
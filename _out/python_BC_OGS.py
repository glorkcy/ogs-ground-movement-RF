# Collection of python boundary condition (BC) classes for OpenGeoSys
import numpy as np

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


import Mesh_settings
ey = Mesh_settings.ey

g = 9.81  #kN/cubic meter


class PorePressure(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, coords, _node_id, _primary_vars):
        if ey==1:
            z = coords[1]
        else:
            z = coords[2]
        
        GWT = -50 + 1 * _t
        distance_GWT_z = GWT - z
        distance_GWT_z = np.asarray(distance_GWT_z)
        depth_under_GWT = np.heaviside([distance_GWT_z], 0) * (distance_GWT_z)
        
        value =  g * depth_under_GWT * 1000 
        return (True, value)

# instantiate the BC objects used by OpenGeoSys
# ---------------------------------------------

PorePressure_py = PorePressure()

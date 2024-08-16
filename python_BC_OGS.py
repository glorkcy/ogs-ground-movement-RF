# Collection of python boundary condition (BC) classes for OpenGeoSys
import numpy as np

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys




g = 9.81  #kN/cubic meter


class PorePressure(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, coords, _node_id, _primary_vars):
        x, y, z = coords
        
        GWT = -50 + 1 * _t
        distance_GWT_y = GWT - y
        distance_GWT_y = np.asarray(distance_GWT_y)
        depth_under_GWT = np.heaviside([distance_GWT_y], 0) * (distance_GWT_y)
        
        value =  g * depth_under_GWT * 1000
        return (True, value)

    

# instantiate the BC objects used by OpenGeoSys
# ---------------------------------------------

PorePressure_py = PorePressure()
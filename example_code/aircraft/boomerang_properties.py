import json
import numpy as np

from atmospheric_functions import statee, gravity_english

class AircraftProperties:
    def __init__(self, V=355., H=24000., cg_shift=[0.0,0.0,0.0], inp_dir='./', **kwargs):

        '''
        Stores and updates aircraft and flight condition specific properties
        
        Parameters
        -----------
        V: float
            airspeed to trim at in ft/s
        H: float
            altitude in feet
        cg_shift: array like
            manual shift in center of gravity location
        inp_dir: string
            directory where aircraft properties input JSON is located
        '''

        # Boomerang properties json input filename
        fn = kwargs.get('filename', 'boomerang_props.json')
        # fn = filename
        prop_dict = json.load(open(inp_dir + fn))
        
        #parse values from JSON input
        self.S_w = prop_dict["geometry"]["S_w"]
        self.b_w = prop_dict["geometry"]["b_w"]
        self.c_w = prop_dict["geometry"]["c_w"]

        # slug-ft^2
        self.Ixx = prop_dict["weight"]["Ixx"]
        self.Iyy = prop_dict["weight"]["Iyy"]
        self.Izz = prop_dict["weight"]["Izz"]
        self.Ixy = prop_dict["weight"]["Ixy"]
        self.Iyx = prop_dict["weight"]["Iyx"]
        self.Ixz = prop_dict["weight"]["Ixz"]
        self.Izx = prop_dict["weight"]["Izx"]
        self.Iyz = prop_dict["weight"]["Iyz"]
        self.Izy = prop_dict["weight"]["Izy"]
        
        self.CG_loc = prop_dict["weight"]["CG"]
        self.cg_shift = cg_shift

        
        self.W = prop_dict["weight"]["W"]
        self.hz = prop_dict["weight"]["hz"]
        self.hy = prop_dict["weight"]["hy"]
        self.hx = 0.0
                
        # some additional flight condition properties
        self.g = gravity_english(H)
        self.mv = self.W/self.g
        _, _, _, _, self.rho, self.a = statee(H)
        _, _, _, _, self.rho_0, self.a_0 = statee(0.0)
        self.nondim_const = 0.5*self.rho*V*V*self.S_w
        self.V = V
        self.H = H
        self.M = self.V/self.a
        
        # shift inertias for manual CG shift
        self.parallel_axis_theorem(s = self.cg_shift, m_cg = self.mv)
        
    def update_propeller_hx(self,hx_sum):
        
        '''
        Updates the propeller angular momentum variable within this class. 
        Values are calculated externally, in the boomerang_aero.py functions.
        '''
        
        self.hx = hx_sum
        
    def parallel_axis_theorem(self, s, m_cg):
        
        '''
        Shifts aircraft moments of inertia to new location using PAT
        
        s: array like
            array defining the center of gravity shift
        m_cg: float
            mass of the vehicle
        '''
        
        I1 = np.array([[self.Ixx, -self.Ixy, -self.Ixz],
                      [-self.Iyx, self.Iyy, -self.Iyz],
                      [-self.Izx, -self.Izy, self.Izz]])

        # 3X3 dentity matrix E
        E = np.identity(3)
        
        # Outer product of s with itself
        sst = np.outer(s, s)
                
        # Add temp to I1
        I2 = I1 + m_cg*(np.dot(s, s)*E - sst)
        
        self.Ixx =  I2[0,0]
        self.Ixy = -I2[0,1]
        self.Ixz = -I2[0,2]
        self.Iyx = -I2[1,0]
        self.Iyy =  I2[1,1]
        self.Iyz = -I2[1,2]
        self.Izx =  -I2[2,0]
        self.Izy =  -I2[2,1]
        self.Izz =  I2[2,2]

if __name__ == "__main__":
    
    test = AircraftProperties()
    
    I1 = np.array([[10.11019833751, 0.0, -1.255109928009],
                   [0.0, 13.05103714647, 0.0],
                   [-1.255109928009, 0.0, 21.89309550680]])
    
    s = -np.array([-0.01607729772902, 0.0, 0.2694770747011])
    
    m = 3.108095017157
    
    I0 = test.parallel_axis_theorem(I1, s, m_cg = m)
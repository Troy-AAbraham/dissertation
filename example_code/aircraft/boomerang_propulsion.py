import numpy as np
import json

class PropModel:
    def __init__(self, inp_dir='./', **kwargs):
        
        '''Pulls coefficients from thrust model JSON and assigns member
        variables.
        
        Values that could be updated in the boomerang_thrust_model.json
        
        dr: float
            propeller diameter in ft
            
        Pp_R and Pp_L: array like
            arrays defining the position (x,y,z) of the propeller component 
            relative to the aircrafts center of gravity in ft
            
        phi,theta,psi: float
            Euler angles defining the orientation of propeller components
            relative to the aircrafts orientation, given in deg in the JSON.
        '''
        
        fn = kwargs.get('fn', 'boomerang_thrust_model.json')
        self.model_coeffs_dict = json.load(open(inp_dir + fn))
        
        # propeller diameter, fixed pitch-to-diameter ratio, moment of inertia
        # locations of right and left propellers
        self.dr = self.model_coeffs_dict["prop_dims"]["dr"]
        self.Kc = self.model_coeffs_dict["prop_dims"]["Kc"]
        self.Ixx = self.model_coeffs_dict["prop_dims"]["Ixx"]
        self.Pp_R = self.model_coeffs_dict["prop_dims"]["Pp_R"]
        self.Pp_L = self.model_coeffs_dict["prop_dims"]["Pp_L"]
        
        # Euler angle defined orientations of propellers relative to aircraft CG
        phi_c = np.deg2rad(self.model_coeffs_dict["prop_dims"]["phi"])
        theta_c = np.deg2rad(self.model_coeffs_dict["prop_dims"]["theta"])
        psi_c = np.deg2rad(self.model_coeffs_dict["prop_dims"]["psi"])
        
        # direction cosine matrix for conversion from propeller to aircraft
        # coordinate system
        self.dir_cos = self.comp_2_body_mat(phi_c, theta_c, psi_c)
        
        # Coefficients defining models to determine the equations for propeller
        # coefficients as a function of pitch-to-diameter ratio. Based on the 
        # model presented in Hunsakers Simulation of Flight
        self.CT_coeffs = self.model_coeffs_dict["CT"]
        self.CP_coeffs = self.model_coeffs_dict["CP"]
        self.CN_a_coeffs = self.model_coeffs_dict["CN_a"]
        self.Cn_a_coeffs = self.model_coeffs_dict["Cn_a"]
        
        self.CT0K0 = self.CT_coeffs["CT0"]["K0"]
        self.CT0K1 = self.CT_coeffs["CT0"]["K1"]
        self.CT0K2 = self.CT_coeffs["CT0"]["K2"]
        
        self.CT1K0 = self.CT_coeffs["CT1"]["K0"]
        self.CT1K1 = self.CT_coeffs["CT1"]["K1"]
        self.CT1K2 = self.CT_coeffs["CT1"]["K2"]

        self.CT2K0 = self.CT_coeffs["CT2"]["K0"]
        self.CT2K1 = self.CT_coeffs["CT2"]["K1"]
        self.CT2K2 = self.CT_coeffs["CT2"]["K2"]
        self.CT2K3 = self.CT_coeffs["CT2"]["K3"]
        
        self.CP0K0 = self.CP_coeffs["CP0"]["K0"]
        self.CP0K1 = self.CP_coeffs["CP0"]["K1"]
        self.CP0K2 = self.CP_coeffs["CP0"]["K2"]
        self.CP0K3 = self.CP_coeffs["CP0"]["K3"]
        self.CP0K4 = self.CP_coeffs["CP0"]["K4"]
        
        self.CP1K0 = self.CP_coeffs["CP1"]["K0"]
        self.CP1K1 = self.CP_coeffs["CP1"]["K1"]
        self.CP1K2 = self.CP_coeffs["CP1"]["K2"]
        self.CP1K3 = self.CP_coeffs["CP1"]["K3"]
        self.CP1K4 = self.CP_coeffs["CP1"]["K4"]
        
        self.CP2K0 = self.CP_coeffs["CP2"]["K0"]
        self.CP2K1 = self.CP_coeffs["CP2"]["K1"]
        self.CP2K2 = self.CP_coeffs["CP2"]["K2"]
        self.CP2K3 = self.CP_coeffs["CP2"]["K3"]
        self.CP2K4 = self.CP_coeffs["CP2"]["K4"]
        
        self.CN1K0 = self.CN_a_coeffs["CN1"]["K0"]
        self.CN1K1 = self.CN_a_coeffs["CN1"]["K1"]
        self.CN1K2 = self.CN_a_coeffs["CN1"]["K2"]
        
        self.CN2K0 = self.CN_a_coeffs["CN2"]["K0"]
        self.CN2K1 = self.CN_a_coeffs["CN2"]["K1"]
        self.CN2K2 = self.CN_a_coeffs["CN2"]["K2"]

        self.CN3K0 = self.CN_a_coeffs["CN3"]["K0"]
        self.CN3K1 = self.CN_a_coeffs["CN3"]["K1"]
        self.CN3K2 = self.CN_a_coeffs["CN3"]["K2"]

        self.Cn1K0 = self.Cn_a_coeffs["Cn1"]["K0"]
        self.Cn1K1 = self.Cn_a_coeffs["Cn1"]["K1"]
        self.Cn1K2 = self.Cn_a_coeffs["Cn1"]["K2"]
        
        self.Cn2K0 = self.Cn_a_coeffs["Cn2"]["K0"]
        self.Cn2K1 = self.Cn_a_coeffs["Cn2"]["K1"]
        self.Cn2K2 = self.Cn_a_coeffs["Cn2"]["K2"]

        self.Cn3K0 = self.Cn_a_coeffs["Cn3"]["K0"]
        self.Cn3K1 = self.Cn_a_coeffs["Cn3"]["K1"]
        self.Cn3K2 = self.Cn_a_coeffs["Cn3"]["K2"]

    def calc_Kc_coeffs(self):
        '''
        Solves for the coefficients that will define the equations for thrust, 
        power, normal force, and yawing moment coefficients. Based on the methods
        outline in Simulation of Flight by Hunsaker (unpublished as of 3/30/25)
        
        '''
        Kc = self.Kc
        Kc2 = Kc*Kc
        Kc3 = Kc2*Kc
        Kc4 = Kc3*Kc

        self.CT0 = self.CT0K2*Kc2 + self.CT0K1*Kc + self.CT0K0
        self.CT1 = self.CT1K2*Kc2 + self.CT1K1*Kc + self.CT1K0
        self.CT2 = self.CT2K3*Kc3 + self.CT2K2*Kc2 + self.CT2K1*Kc + self.CT2K0
        
        self.CP0 = self.CP0K4*Kc4 + self.CP0K3*Kc3 + self.CP0K2*Kc2 + self.CP0K1*Kc + self.CP0K0
        self.CP1 = self.CP1K4*Kc4 + self.CP1K3*Kc3 + self.CP1K2*Kc2 + self.CP1K1*Kc + self.CP1K0
        self.CP2 = self.CP2K4*Kc4 + self.CP2K3*Kc3 + self.CP2K2*Kc2 + self.CP2K1*Kc + self.CP2K0
        
        self.CN1 = self.CN1K2*Kc2 + self.CN1K1*Kc + self.CN1K0
        self.CN2 = self.CN2K2*Kc2 + self.CN2K1*Kc + self.CN2K0
        self.CN3 = self.CN3K2*Kc2 + self.CN3K1*Kc + self.CN3K0
        
        self.Cn1 = self.Cn1K2*Kc2 + self.Cn1K1*Kc + self.Cn1K0
        self.Cn2 = self.Cn2K2*Kc2 + self.Cn2K1*Kc + self.Cn2K0
        self.Cn3 = self.Cn3K2*Kc2 + self.Cn3K1*Kc + self.Cn3K0
    
    
    def calc_coeffs(self,J):
        '''
        Solves for the thrust, power, normal force derivative, and yawing moment
        derivative as a function of advance ratio
        
        Parameters
        -----------
        
        J: float
            Propeller advance ratio - dimensionless
            
        Returns
        -----------
        
        CT: float
            Thrust coefficient
            
        CP: float
            Power coefficient
            
        CN_a: float
            Normal force coefficient derivative
            
        Cn_a: float
            Yawing moment coefficient derivative
        '''
        
        CT = self.CT0 + self.CT1*J + self.CT2*J*J
        CP = self.CP0 + self.CP1*J + self.CP2*J*J
        CN_a = self.CN1*J + self.CN2*J*J + self.CN3*J*J*J
        Cn_a = self.Cn1*J + self.Cn2*J*J + self.Cn3*J*J*J
        
        return CT,CP,CN_a,Cn_a
    
    def comp_2_body_mat(self,phi_c,theta_c,psi_c):
        '''
        Creates the direction cosine matrix for the conversion between the 
        propeller coordinate system and the aircraft coordinate system.
        
        Uses the form that is propeller to aircraft conversion.
        
        Parameters
        -----------
        phi_c,theta_c,psi_c: floats
            Euler angles [rad]
            
        Returns
        -----------
        dir_cosine: array-like (3x3 matrix)
        
        '''

        SP,ST,SPS = np.sin(phi_c), np.sin(theta_c), np.sin(psi_c)
        CP,CT,CPS = np.cos(phi_c), np.cos(theta_c), np.cos(psi_c)
        
        dir_cosine = np.array([[CT*CPS, SP*ST*CPS-CP*SPS, CP*ST*CPS + SP*SPS],
                             [CT*SPS, SP*ST*SPS + CP*CPS, CP*ST*SPS - SP*CPS],
                             [-ST, SP*CT, CP*CT]])
        return dir_cosine
    
    def prop_velocity(self,V_body,omega_body,Pp_body):
        '''
        Calculate the local velocity at the propeller as a result of aircraft
        translational velocity, rotational velocity, and propeller position.
        
        Parameters
        -----------
        V_body: array_like
            vehicle body-fixed velocity vector
            
        omega_body: array_like
            vehicle body-fixed roational velocity vector
            
        Pp_body: array_like
            body-fixed position of the component relative to vehicle CG
            
        Returns
        -----------
        V_comp: array-like
            Velocity vector in the propellers coordinate system
            
        V_mag: float
            Velocity magnitude at the propeller
            
        u_comp: array-like
            Unit velocity vector in the propellers coordinate system
            
        '''
        
        V_comp = np.matmul(np.linalg.inv(self.dir_cos),(V_body - np.cross(Pp_body,omega_body)))
        V_mag = np.sqrt(V_comp[0]*V_comp[0] + V_comp[1]*V_comp[1] + V_comp[2]*V_comp[2])
        u_comp = V_comp/V_mag
        
        return V_comp, V_mag, u_comp
    
    def solve_J(self,T,Vinf,rho):
        '''
        Find the advance ratio (J) that will produce the desired thrust for the
        given airspeed and air density.
        
        Parameters
        -----------
        
        T: float
            desired thrust
            
        Vinf: float
            airspeed
            
        rho: float
            airdensity
            
        Returns
        -----------
        
        J: float
            Propeller advance ratio - dimensionless
        '''
        # be aware of the sign of this solution
        J = (-self.CT1 - np.sqrt(self.CT1*self.CT1 - 4*(self.CT2 - (T/(rho*Vinf*Vinf*self.dr*self.dr)))*self.CT0))/(2*(self.CT2 - (T/(rho*Vinf*Vinf*self.dr*self.dr))))
        
        return J 

    def calc_FM(self,V_body,omega_body,Pp_body,rho,omega,delta=1):
        '''
        Finds the velocity components at the propeller. Then finds the 
        propeller force and moment coefficients as a function of air density, 
        propeller rotation rate, and alpha. Then dimensionalizes those coefficients.
        
        Parameters
        -----------
        
        V_body: array_like
            vehicle body-fixed velocity vector
            
        omega_body: array_like
            vehicle body-fixed roational velocity vector
            
        Pp_body: array_like
            body-fixed position of the component relative to vehicle CG
            
        rho: float
            atmospheric density
            
        omega: float
            propeller rotation rate, driven by throttle setting tau
            
        delta: int (1 or -1)
            defines rotor as right handed (1) or left handed (-1)
            
        Returns
        -----------
        F_vec: array_like
            vehicle body-fixed force vector
            
        M_vec: array_like
            vehicle body-fixed moment vector
            
        '''
        
        # find velocity in the propeller coordinate system
        V_comp, V_mag_comp, u_comp = self.prop_velocity(V_body,omega_body,Pp_body)
        Vinf = V_mag_comp #
        # propeller angle of attack
        alpha = np.arccos(u_comp[0])
        
        #calculate advance ratio
        J = (2*np.pi*Vinf)/(omega*self.dr)
        
        # propeller force and moment coefficients as a function of advance ratio
        CT = self.CT0 + self.CT1*J + self.CT2*J*J
        CP = self.CP0 + self.CP1*J + self.CP2*J*J
        CN_a = self.CN1*J + self.CN2*J*J + self.CN3*J*J*J
        Cn_a = self.Cn1*J + self.Cn2*J*J + self.Cn3*J*J*J
        
        # dimensionalize forces and moments
        T = rho*((omega/(2*np.pi))**2)*(self.dr**4)*CT
        P = rho*((omega/(2*np.pi))**3)*(self.dr**5)*CP
        N = rho*((omega/(2*np.pi))**2)*(self.dr**4)*CN_a*alpha
        n = rho*((omega/(2*np.pi))**2)*(self.dr**5)*Cn_a*alpha
        ell = P/omega
        
        # determine normal force vector in propeller coordinate system
        if V_comp[1] == 0.0 and V_comp[2] == 0.0:
            U_N = np.array([0.0,0.0,0.0])
        else:
            U_N = np.array([0.0, -V_comp[1], -V_comp[2]])/np.sqrt(V_comp[1]*V_comp[1] + V_comp[2]*V_comp[2])
        
        # sum forces in propeller coordinate system
        F_vec = np.array([T,0.0,0.0]) + N*U_N
        # convert back to aircraft coordinate system
        F_vec = np.matmul(self.dir_cos,F_vec)
        
        # sum moments in propeller coordinate system
        M_vec = np.array([-delta*ell, 0.0, 0.0]) - delta*n*U_N
        # convert back to aircraft coordinate system
        M_vec = np.matmul(self.dir_cos,M_vec)

        return F_vec, M_vec
    
    def get_thrust(self,tau,V,V_vec,omega_vec,rho,cg_shift):
        '''
        Generates total forces and momnets from both propellers.tau will control
        the rotation rate (omega_R) of the right engine. The
        rotation rate of the left engine is then determined in order to maintain
        no yawing moment resulting from the thrust of both engines.
        
        Parameters
        -----------
        
        tau: float
            throttle setting, currently not limited
            
        V: float
            vehicle airspeed
            
        V_vec: array_like
            vehicle body-fixed velocity vector
            
        omega_body: array_like
            vehicle body-fixed roational velocity vector
            
        rho: float
            atmospheric density
            
        Returns
        -----------
        F_vec: array_like
            Vehicle body-fixed propeller force vector
            
        M_vec: array_like
            vehicle body-fixed propeller moment vector
            
        '''
        # positions of the two propellers after a cg shift in the body coordinate system
        Pp_R = [self.Pp_R[0] - cg_shift[0], self.Pp_R[1] - cg_shift[1], self.Pp_R[2] - cg_shift[2]]
        Pp_L =  [self.Pp_L[0] - cg_shift[0], self.Pp_L[1] - cg_shift[1], self.Pp_L[2] - cg_shift[2]]
        
        # max and min propeller rotation rates, chosen somewhat arbitrarily
        # in RPMs
        omega_min = 2000
        omega_max = 3000
        
        # throttle setting tau controls right propeller rotation rate
        omega_R = omega_min + (omega_max - omega_min)*tau
        # convert to rad/s
        omega_R = omega_R*2*np.pi/60.

        # calculate forces and moments on right propeller
        F_vec_R, M_vec_R = self.calc_FM(V_body=V_vec,omega_body=omega_vec,Pp_body=Pp_R,rho=rho,omega=omega_R)
        # moments produced by forces acting away from CG
        M_vec_R_pos = np.cross(Pp_R,F_vec_R)
        # sum total moments to aircraft
        M_vec_R = M_vec_R + M_vec_R_pos
                
        omega_L = omega_R
        
        # should produce T_L thrust in addition to the off axis forces
        F_vec_L, M_vec_L = self.calc_FM(V_body=V_vec,omega_body=omega_vec,Pp_body=Pp_L,rho=rho,omega=omega_L)
        M_vec_L_pos = np.cross(Pp_L,F_vec_L)
        
        M_vec_L = M_vec_L + M_vec_L_pos
        
        # sum right and left propeller moments
        F_vec = F_vec_R + F_vec_L
        M_vec = M_vec_R + M_vec_L 
        
        # propeller angular momentum 
        hx_R = omega_R*self.Ixx
        hx_L = omega_L*self.Ixx
        
        # rotate to aircraft coordinates
        h_R = np.matmul(self.dir_cos,np.array([hx_R, 0.0, 0.0]))
        h_L = np.matmul(self.dir_cos,np.array([hx_L, 0.0, 0.0]))
        
        return F_vec, M_vec,hx_R,hx_L
      

if __name__ == "__main__":
    from aero_component_helper import solve_Vb_vector
    
    tau = 0.7
    V = 355 #ft/s
    alpha = np.deg2rad(5.0)
    beta = np.deg2rad(0.0)
    p = np.deg2rad(0.0)
    q = np.deg2rad(0.0)
    r = np.deg2rad(0.0)
    
    cg_shift = [0.0,0.0,0.0]
    
    rho = 0.0023769
    
    V_vec = solve_Vb_vector(alpha,beta,V)
    omega_vec = np.array([p,q,r])
    
    prop_class = PropModel()
    prop_class.calc_Kc_coeffs()
    
    
    F_thrust_vec, M_thrust_vec, hx_R, hx_L = prop_class.get_thrust(tau,V,V_vec,omega_vec,rho, cg_shift)
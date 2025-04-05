import sys
import numpy as np
from scipy.linalg import eig
import matplotlib.pyplot as plt

aero_directory = '../aircraft/'
sys.path.insert(1, aero_directory)

from boomerang_aero import BoomerangAero
from boomerang_properties import AircraftProperties
from trim_functions import solve_trim
from dynamic_derivatives import solveDerivatives

class dynamicAnalysis:
            
    def __init__(self,  write_output = False, output_filename = 'dynamic_output.txt',
                 shss=False, compressible=False, stall=False, coords_approx = False, derivs_approx = False,
                 cg_shift=[0.0, 0.0, 0.0]):
        '''
        Load all the constant value aircraft parameters and analysis inputs
        
        Parameters
        -----------
        
        write_output: boolean
            flag to write eigenvalues to file
            
        output_filename: string
            eigenvalue output filename
            
        shss: boolean
            sets trim type to steady-heading sideslip, otherwise defaults to SCT
            
        compressible: bool
            unused
            
        stall: bool
            unused
            
        coords_approx: bool
            apply the coordinate system approximation to derivatives and inertia components
            
        derivs_approx: bool
            apply the symmetric derivatve approximation
            
        cg_shift: array-like
               vector defining shift in CG
        '''
        
        self.write_output = write_output
        self.output_filename = output_filename
        
        self.cg_shift = cg_shift
        self.shss = shss
        self.compressible = compressible # not currently used
        self.stall = stall # not currently used
        self.coords_approx = coords_approx
        self.derivs_approx = derivs_approx
        
        self.Gamma = 0.5 # relaxation factor for the trim algorithm
        
        # initialize the aerodynamic model
        self.aeroModel = BoomerangAero(inp_dir=aero_directory)
            
    def update_aircraft_properties(self, V, H):
        
        '''
        Updates aircraft property object and internal values
        
        Parameters
        -----------
        
        V: float
            total airspeed, ft/s
            
        H: float
            altitude, ft
        '''
        
        self.aircraft_properties = AircraftProperties(V = V, H = H, cg_shift = self.cg_shift, Gamma = self.Gamma, inp_dir = aero_directory)
        
        # update operating properties
        self.V = self.aircraft_properties.V
        self.H =self.aircraft_properties.H
        self.g = self.aircraft_properties.g
        self.nondim_const = self.aircraft_properties.nondim_const
        self.rho = self.aircraft_properties.rho
        self.rho_0 = self.aircraft_properties.rho_0
        self.M = self.aircraft_properties.M
        self.a = self.aircraft_properties.a
        self.a_0 = self.aircraft_properties.a_0
        
        # update aircraft geometric properties
        self.bw = self.aircraft_properties.b_w
        self.cw = self.aircraft_properties.c_w
        self.Sw = self.aircraft_properties.S_w
        # update aircraft mass and inertia properties
        self.W = self.aircraft_properties.W
        
        self.Ixxb = self.aircraft_properties.Ixx
        self.Iyyb = self.aircraft_properties.Iyy
        self.Izzb = self.aircraft_properties.Izz
        self.Ixyb = self.aircraft_properties.Ixy
        self.Ixzb = self.aircraft_properties.Ixz
        self.Iyzb = self.aircraft_properties.Iyz

        # self.aircraft_properties.hx = 0.0 # for testing
        self.hxb = self.aircraft_properties.hx
        self.hyb = self.aircraft_properties.hy
        self.hzb = self.aircraft_properties.hz
        
    def solve_equilibrium_state(self, V, H, gamma, phi):
        
        '''
        Given airspeed, altitude, CG shift, climb angle,
        and bank angle, solve equilibrium trim condition.
        
        Parameters
        -----------
        
        V: float
            total airspeed, ft/s
            
        H: float
            altitude, ft
            
        gamma: float
            climb angle, rad
            
        phi: float
            bank angle, rad
            
        '''
        # update aircraft properties
        self.update_aircraft_properties(V, H)
        
        #solve for trim solution
        self.solution = solve_trim(aero_model = self.aeroModel, aircraft_props=self.aircraft_properties,
                              gamma = gamma, phi = phi, cg_shift=self.cg_shift,
                              shss = self.shss, compressible = self.compressible, stall = self.stall)
        
        #update after trim solution to get correct angular momentum from propellers
        self.hxb = self.aircraft_properties.hx
        
        #parse trim solution values
        tau, alpha, beta, da, de, dr = self.solution.x
        u, v, w, p, q, r, phi, theta = self.solution.states
        
        FX, FY, FZ, Mx, My, Mz = self.solution.FM_dim
        [CL, CS, CD, Cl, Cm, Cn] = self.solution.FM
        
        # set equilibrium velocity vector depending on approximations used
        # and trim solution
        if self.coords_approx:
            self.eq_velo = np.array([self.V,0.0,0.0])
        else:
            self.eq_velo = np.array([u,v,w])
        
        # set equilibrium rotation rate, euler angles, control inputs, and FM vectors
        self.eq_rot = np.array([p,q,r])
        self.eq_euler = np.array([phi,theta])
        self.eq_inputs = np.array([tau, alpha, beta, da, de, dr])
        self.eq_FM_wind = np.array([CL, CS, CD, Cl, Cm, Cn])
        self.eq_FM = np.array([FX, FY, FZ, Mx, My, Mz])
        
        # set equilibrium aerodynamic angles
        self.alpha = alpha
        self.beta = beta
    
    def solve_derivatives(self, print_latex_tables = False):
        
        '''        
        Solves for the required body-fixed derivatives.
        
        Parameters
        -----------
        
        print_latex_tables: boolean
            flag to print derivative solutions
            
        
        '''
                
        # initialize derivative class
        derivs = solveDerivatives(aeroModel = self.aeroModel, aircraft_properties = self.aircraft_properties, cg_shift = self.cg_shift,
                    trim_solution = self.solution, compressible = self.compressible, stall = self.stall,
                    coords_approx = self.coords_approx, derivs_approx = self.derivs_approx)
        
        # solve derivative solution
        deriv_solution = derivs.solve_derivs()
        # update local derivative solution
        self.set_deriv_solution(deriv_solution)
        # print derivatives
        derivs.print_derivs(print_control = True, print_latex_tables = print_latex_tables)
            
    def set_deriv_solution(self,deriv_solution):
        
        '''Sets derivatives in the current class.
        
        Parameters
        -----------
        
        deriv_solution: object
            parse derivatives from the deriv_solution object
            
        '''
        
        self.Fxb_u = deriv_solution.Fxb_u
        self.Fxb_v = deriv_solution.Fxb_v
        self.Fxb_w = deriv_solution.Fxb_w
        self.Fxb_p = deriv_solution.Fxb_p
        self.Fxb_q = deriv_solution.Fxb_q
        self.Fxb_r = deriv_solution.Fxb_r
        
        self.Fyb_u = deriv_solution.Fyb_u
        self.Fyb_v = deriv_solution.Fyb_v
        self.Fyb_w = deriv_solution.Fyb_w
        self.Fyb_p = deriv_solution.Fyb_p
        self.Fyb_q = deriv_solution.Fyb_q
        self.Fyb_r = deriv_solution.Fyb_r
        
        self.Fzb_u = deriv_solution.Fzb_u
        self.Fzb_v = deriv_solution.Fzb_v
        self.Fzb_w = deriv_solution.Fzb_w
        self.Fzb_p = deriv_solution.Fzb_p
        self.Fzb_q = deriv_solution.Fzb_q
        self.Fzb_r = deriv_solution.Fzb_r
                
        self.Mxb_u = deriv_solution.Mxb_u
        self.Mxb_v = deriv_solution.Mxb_v 
        self.Mxb_w = deriv_solution.Mxb_w
        self.Mxb_p = deriv_solution.Mxb_p
        self.Mxb_q = deriv_solution.Mxb_q
        self.Mxb_r = deriv_solution.Mxb_r 
        
        self.Myb_u = deriv_solution.Myb_u
        self.Myb_v = deriv_solution.Myb_v
        self.Myb_w = deriv_solution.Myb_w
        self.Myb_p = deriv_solution.Myb_p
        self.Myb_q = deriv_solution.Myb_q
        self.Myb_r = deriv_solution.Myb_r
        
        self.Mzb_u = deriv_solution.Mzb_u
        self.Mzb_v = deriv_solution.Mzb_v
        self.Mzb_w = deriv_solution.Mzb_w
        self.Mzb_p = deriv_solution.Mzb_p
        self.Mzb_q = deriv_solution.Mzb_q
        self.Mzb_r = deriv_solution.Mzb_r
        
        self.Fzb_wdot = deriv_solution.Fzb_wdot
        self.Myb_wdot = deriv_solution.Myb_wdot

        self.Fxb_udot = deriv_solution.Fxb_udot
        self.Fxb_vdot = deriv_solution.Fxb_vdot
        self.Fxb_wdot = deriv_solution.Fxb_wdot
        self.Fyb_udot = deriv_solution.Fyb_udot
        self.Fyb_vdot = deriv_solution.Fyb_vdot
        self.Fyb_wdot = deriv_solution.Fyb_wdot
        self.Fzb_vdot = deriv_solution.Fzb_vdot
        
        self.Fzb_udot = deriv_solution.Fzb_udot
        self.Myb_udot = deriv_solution.Myb_udot
        
        self.Fxb_da = deriv_solution.Fxb_da
        self.Fyb_da = deriv_solution.Fyb_da
        self.Fzb_da = deriv_solution.Fzb_da
        self.Mxb_da = deriv_solution.Mxb_da
        self.Myb_da = deriv_solution.Myb_da
        self.Mzb_da = deriv_solution.Mzb_da
        
        self.Fxb_de = deriv_solution.Fxb_de
        self.Fyb_de = deriv_solution.Fyb_de
        self.Fzb_de = deriv_solution.Fzb_de
        self.Mxb_de = deriv_solution.Mxb_de
        self.Myb_de = deriv_solution.Myb_de
        self.Mzb_de = deriv_solution.Mzb_de
        
        self.Fxb_dr = deriv_solution.Fxb_dr
        self.Fyb_dr = deriv_solution.Fyb_dr
        self.Fzb_dr = deriv_solution.Fzb_dr
        self.Mxb_dr = deriv_solution.Mxb_dr
        self.Myb_dr = deriv_solution.Myb_dr
        self.Mzb_dr = deriv_solution.Mzb_dr
        
        self.Fxb_tau = deriv_solution.Fxb_tau
        self.Fyb_tau = deriv_solution.Fyb_tau
        self.Fzb_tau = deriv_solution.Fzb_tau
        self.Mxb_tau = deriv_solution.Mxb_tau
        self.Myb_tau = deriv_solution.Myb_tau
        self.Mzb_tau = deriv_solution.Mzb_tau
        
        self.Ixxb = deriv_solution.Ixxb
        self.Iyyb = deriv_solution.Iyyb
        self.Izzb = deriv_solution.Izzb
        self.Ixyb = deriv_solution.Ixyb
        self.Ixzb = deriv_solution.Ixzb
        self.Iyzb = deriv_solution.Iyzb
    
    def generate_control_matrix(self):
        
        '''Generates linear system control matrix'''

        self.Control_matrix = np.array([[self.Fxb_da,                 self.Fxb_de,                 self.Fxb_dr,                 self.Fxb_tau],
                                        [self.Fyb_da,                 self.Fyb_de,                 self.Fyb_dr,                 self.Fyb_tau],
                                        [self.Fzb_da,                 self.Fzb_de,                 self.Fzb_dr,                 self.Fzb_tau],
                                        [self.Mxb_da,                 self.Mxb_de,                 self.Mxb_dr,                 self.Mxb_tau],
                                        [self.Myb_da,                 self.Myb_de,                 self.Myb_dr,                 self.Myb_tau],
                                        [self.Mzb_da,                 self.Mzb_de,                 self.Mzb_dr,                 self.Mzb_tau],
                                        [0.0,                         0.0,                         0.0,                         0.0],
                                        [0.0,                         0.0,                         0.0,                         0.0],
                                        [0.0,                         0.0,                         0.0,                         0.0],
                                        [0.0,                         0.0,                         0.0,                         0.0],
                                        [0.0,                         0.0,                         0.0,                         0.0],
                                        [0.0,                         0.0,                         0.0,                         0.0]])
    
    def solve_dynamics_system(self, print_results = True, remove_xyz = False, norm_types = False):

        '''
        Solves the LTI system and report eigensolution. Uses the dimensional
        linearized equations of motion presented in my dissertation.
        
        Parameters
        -----------
        
        print_results: boolean
            flag for printing eigenvalue and eigenvector solutions
            
        remove_xyz: boolean
            flag for removing the xyz components from the LTI
            
        norm_types: boolean
            flag for normalizing the eigenvector components relative to the seperate
            state variable types (velocities, rotation rates, positions, orientations)

        '''
        
        # set asymmetrix producats of inertia to zero if symmetric aircraft assumption is applied
        if self.derivs_approx == True:
            self.Ixyb = 0.0
            self.Iyzb = 0.0
        
        # set equilibrium body-fixed velocities
        u_o = self.eq_velo[0]
        v_o = self.eq_velo[1]
        w_o = self.eq_velo[2]
        
        # set equilibrium bank and elevation angles
        phi_o = (self.eq_euler[0])
        theta_o = (self.eq_euler[1])
        
        '''INTERNAL CALCULATIONS'''
        # solve for the trig values
        ST_o = np.sin(theta_o)
        CT_o = np.cos(theta_o)
        TT_o = np.tan(theta_o)
        
        SP_o = np.sin(phi_o)
        CP_o = np.cos(phi_o)
        
        if self.shss!= True:
            '''If not SHSS, used SCT equations to get pqr, these should match the trim solutions.
            Using this form allows the coordinate system assumption to influence these values.'''
            p_o = (-self.g*SP_o*ST_o*CT_o)/(w_o*ST_o + u_o*CP_o*CT_o)
            q_o = (self.g*SP_o*SP_o*CT_o*CT_o)/(w_o*ST_o + u_o*CP_o*CT_o)
            r_o = (self.g*SP_o*CP_o*CT_o*CT_o)/(w_o*ST_o + u_o*CP_o*CT_o)
        else:
            '''If SHSS, set equal to trim solution values. These should always be zero.
            They are not zero if you used the equations above in a SHSS.'''
            p_o = self.eq_rot[0]
            q_o = self.eq_rot[1]
            r_o = self.eq_rot[2]
        
        self.p0 = p_o
        self.q0 = q_o
        self.r0 = r_o

        print('\n')
        
        # solve aircraft mass, represented as ratio in the linear EOM
        W_g = self.W/self.g
        
        # set B Matrix
        B_matrix = np.array([[W_g - self.Fxb_udot, 0.0,      -self.Fxb_wdot,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, W_g,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [     -self.Fzb_udot, 0.0, W_g - self.Fzb_wdot,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,  self.Ixxb, -self.Ixyb, -self.Ixzb, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [     -self.Myb_udot, 0.0,      -self.Myb_wdot, -self.Ixyb,  self.Iyyb, -self.Iyzb, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0, -self.Ixzb, -self.Iyzb,  self.Izzb, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
        
        AMxp = self.Mxb_p + self.Ixzb*q_o - self.Ixyb*r_o
        AMxq = self.Mxb_q - self.hzb + (self.Iyyb - self.Izzb)*r_o + 2*self.Iyzb*q_o + self.Ixzb*p_o
        AMxr = self.Mxb_r + self.hyb + (self.Iyyb - self.Izzb)*q_o - 2*self.Iyzb*r_o - self.Ixyb*p_o
        AMyp = self.Myb_p + self.hzb + (self.Izzb - self.Ixxb)*r_o - 2*self.Ixzb*p_o - self.Iyzb*q_o
        AMyq = self.Myb_q + self.Ixyb*r_o - self.Iyzb*p_o
        AMyr = self.Myb_r - self.hxb + (self.Izzb - self.Ixxb)*p_o + 2*self.Ixzb*r_o + self.Ixyb*q_o
        AMzp = self.Mzb_p - self.hyb + (self.Ixxb - self.Iyyb)*q_o + 2*self.Ixyb*p_o + self.Iyzb*r_o
        AMzq = self.Mzb_q + self.hxb + (self.Ixxb - self.Iyyb)*p_o - 2*self.Ixyb*q_o - self.Ixzb*r_o
        AMzr = self.Mzb_r + self.Iyzb*p_o - self.Ixzb*q_o

        AxP = v_o*(CP_o*ST_o) + w_o*(-SP_o*ST_o)
        AxT = -u_o*ST_o + v_o*SP_o*CT_o+ w_o*CP_o*CT_o
        AxO = - v_o*(CP_o) + w_o*(SP_o)
        AyP = v_o*(- SP_o) - w_o*(CP_o)
        AyT = 0.0
        AyO = u_o*CT_o + v_o*(SP_o*ST_o) + w_o*(CP_o*ST_o)
        AzP = v_o*CP_o*CT_o - w_o*SP_o*CT_o
        AzT = -u_o*CT_o - v_o*SP_o*ST_o - w_o*CP_o*ST_o
        AzO = 0.0
        
        A441 = 0
        A442 = (q_o*SP_o + r_o*CP_o)*(1/(CT_o*CT_o))
        A443 = (-q_o*SP_o - r_o*CP_o)
        A444 = 0
        A445 = (q_o*SP_o + r_o*CP_o)*(TT_o/CT_o)
        
        # set A Matrix
        A_matrix = np.array([[          self.Fxb_u,       self.Fxb_v + W_g*r_o,       self.Fxb_w - W_g*q_o,           self.Fxb_p, self.Fxb_q - W_g*w_o, self.Fxb_r + W_g*v_o,       0.0, 0.0, 0.0,               0.0,      -self.W*CT_o, 0.0],
                             [self.Fyb_u - W_g*r_o,                 self.Fyb_v,       self.Fyb_w + W_g*p_o, self.Fyb_p + W_g*w_o,           self.Fyb_q, self.Fyb_r - W_g*u_o,       0.0, 0.0, 0.0,  self.W*CP_o*CT_o, -self.W*SP_o*ST_o, 0.0],
                             [self.Fzb_u + W_g*q_o,       self.Fzb_v - W_g*p_o,                 self.Fzb_w, self.Fzb_p - W_g*v_o, self.Fzb_q + W_g*u_o,           self.Fzb_r,       0.0, 0.0, 0.0, -self.W*SP_o*CT_o, -self.W*CP_o*ST_o, 0.0],
                             [          self.Mxb_u,                 self.Mxb_v,                 self.Mxb_w,                 AMxp,                 AMxq,                 AMxr,       0.0, 0.0, 0.0,               0.0,               0.0, 0.0],
                             [          self.Myb_u,                 self.Myb_v,                 self.Myb_w,                 AMyp,                 AMyq,                 AMyr,       0.0, 0.0, 0.0,               0.0,               0.0, 0.0],
                             [          self.Mzb_u,                 self.Mzb_v,                 self.Mzb_w,                 AMzp,                 AMzq,                 AMzr,       0.0, 0.0, 0.0,               0.0,               0.0, 0.0],
                             [                CT_o,                  SP_o*ST_o,                  CP_o*ST_o,                  0.0,                  0.0,                  0.0,       0.0, 0.0, 0.0,               AxP,               AxT, AxO],
                             [                   0,                       CP_o,                      -SP_o,                  0.0,                  0.0,                  0.0,       0.0, 0.0, 0.0,               AyP,               AyT, AyO],
                             [               -ST_o,                  SP_o*CT_o,                  CP_o*CT_o,                  0.0,                  0.0,                  0.0,       0.0, 0.0, 0.0,               AzP,               AzT, AzO],
                             [                 0.0,                        0.0,                        0.0,                  1.0,            SP_o*TT_o,            CP_o*TT_o,       0.0, 0.0, 0.0,              A441,              A442, 0.0],
                             [                 0.0,                        0.0,                        0.0,                  0.0,                 CP_o,                -SP_o,       0.0, 0.0, 0.0,              A443,               0.0, 0.0],
                             [                 0.0,                        0.0,                        0.0,                  0.0,            SP_o/CT_o,            CP_o/CT_o,       0.0, 0.0, 0.0,              A444,              A445, 0.0]])
        
        
        # remove xyz components from A and B matrices
        if remove_xyz == True:
            A_matrix = np.delete(A_matrix, [6,7,8], axis=0)
            A_matrix = np.delete(A_matrix, [6,7,8], axis=1)
            
            B_matrix = np.delete(B_matrix, [6,7,8], axis=0)
            B_matrix = np.delete(B_matrix, [6,7,8], axis=1)
                    
        # calculate C matrix
        C_matrix = np.matmul(np.linalg.inv(B_matrix),A_matrix)
        self.B_matrix = B_matrix
        self.A_matrix = A_matrix
        self.C_matrix = C_matrix
        
        # solve for eigenvalues and eigenvectors
        self.eigvals, self.eigvecs = eig(C_matrix)
        self.eigvecs_old = np.copy(self.eigvecs)
        
        # print A B and C matrices
        print('A-Matrix')
        print('\n'.join([''.join(['{:>16.6f}'.format(item) for item in row]) 
                for row in A_matrix]))
        
        print('\n')
        print('B-Matrix')
        print('\n'.join([''.join(['{:>16.6f}'.format(item) for item in row]) 
                for row in B_matrix]))
        
        print('\n')
        print('C-Matrix')
        print('\n'.join([''.join(['{:>16.6f}'.format(item) for item in row]) 
                for row in C_matrix]))    
        
        
        #normalize eigenvectors relative to the largest in each array

        for i in range(len(self.eigvecs[0,:])):
            
            index_max = np.argmax(np.abs(self.eigvecs[:,i]))
            
            '''MULTIPLYING BY CC JUST SHIFTS THE PHASE ANGLES TO BE RELATIVE TO THAT COMPONENT'''
            # index_max = 8 
            
            cc = np.conj(self.eigvecs[index_max,i])
            
            new_vec = cc*self.eigvecs[:,i]
            
            new_vec = new_vec / np.sqrt(np.sum(np.square(np.abs(new_vec))))
            
            self.eigvecs[:,i] = new_vec

        
        # normalize the eigenvectors within their specific state variable types
        if norm_types == True:
            
            if remove_xyz == True:
                
                for i in range(len(self.eigvecs[:,0])):
                    # velocities
                    self.eigvecs[0:3,i] = self.eigvecs[0:3,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[0:3,i]))))
                    # rotatation rates
                    self.eigvecs[3:6,i] = self.eigvecs[3:6,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[3:6,i]))))
                    # orientation
                    self.eigvecs[6:9,i] = self.eigvecs[6:9,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[6:9,i]))))
            else:
                
                for i in range(len(self.eigvecs[:,0])):
                    # velocities
                    self.eigvecs[0:3,i] = self.eigvecs[0:3,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[0:3,i]))))
                    # rotation rates
                    self.eigvecs[3:6,i] = self.eigvecs[3:6,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[3:6,i]))))
                    # earth fixed positions
                    self.eigvecs[6:9,i] = self.eigvecs[6:9,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[6:9,i]))))
                    # orientation
                    self.eigvecs[9:12,i] = self.eigvecs[9:12,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[9:12,i]))))
        
        # sort eigenvalues according to type first (real and imaginary) and then magnitude
        self.eigvals, i_sort = self.sort_eigenvalues_and_indices(self.eigvals)
        # match eigenvectors to sorted eigenvalues
        self.eigvecs = self.eigvecs[:,i_sort]
        # calculate eigenvector amplitudes and phase
        self.amps = np.abs(self.eigvecs)
        self.phase = np.rad2deg(np.arctan2(np.imag(self.eigvecs),np.real(self.eigvecs)))
        
        # parse and generate dynamic stability parameters
        self.eigreal = np.real(self.eigvals[:])
        self.eigimag = np.imag(self.eigvals[:])
        self.sigma = -np.real(self.eigvals[:])
        self.omegad = np.abs(np.imag(self.eigvals))
        self.period = 2.0*np.pi/self.omegad
        
        # print results to terminal
        if print_results == True:
            # print eigenvalues
            print('\nEigenvalues') 
            print('\n'.join('{:>32.12f}'.format(item) for item in self.eigvals))
                    
            if self.write_output != False:
                with open(self.output_filename, 'w') as export_handle:
                    for i in range(len(self.eigvals)):
                        print('{:<16.12f}{:<16.12f}'.format(np.real(self.eigvals[i]), np.imag(self.eigvals[i])), file=export_handle)       
            
            # print eigenvalues, related oscillatory properties, and eigenvectors
            for i in range(len(self.eigvecs[0,:])):
                
                print('\n')
                print('{:>24}'.format('Eigenvalue'), '{:<18.12f}'.format(self.eigvals[i]))
                print('{:>24}'.format('Real'), '{:<18.12f}'.format(self.eigreal[i]))
                print('{:>24}'.format('Imaginary'), '{:<18.12f}'.format(self.eigimag[i]))
                print('{:>24}'.format('Period'), '{:<18.12f}'.format(self.period[i]))
                print('{:>24}'.format('Damping'), '{:<18.12f}'.format(self.sigma[i]))
                
                if len(self.eigvecs[:,0]) == 12:
                    # print("\nEigenvectors:")
                    print('{:>28}'.format('component'),  '{:>28}'.format('eigenvector'), '{:>28}'.format('phase'), '{:>28}'.format('amplitude'))
                    print('{:>28}'.format('\u0394u'),  '{:>28.10f}'.format(self.eigvecs[0,i]), '{:>28.12f}'.format(self.phase[0,i]), '{:>28.12f}'.format(self.amps[0,i]))
                    print('{:>28}'.format('\u0394v'),  '{:>28.10f}'.format(self.eigvecs[1,i]), '{:>28.12f}'.format(self.phase[1,i]), '{:>28.12f}'.format(self.amps[1,i]))
                    print('{:>28}'.format('\u0394w'),  '{:>28.10f}'.format(self.eigvecs[2,i]), '{:>28.12f}'.format(self.phase[2,i]), '{:>28.12f}'.format(self.amps[2,i]))
                    print('{:>28}'.format('\u0394p'),  '{:>28.10f}'.format(self.eigvecs[3,i]), '{:>28.12f}'.format(self.phase[3,i]), '{:>28.12f}'.format(self.amps[3,i]))
                    print('{:>28}'.format('\u0394q'),  '{:>28.10f}'.format(self.eigvecs[4,i]), '{:>28.12f}'.format(self.phase[4,i]), '{:>28.12f}'.format(self.amps[4,i]))
                    print('{:>28}'.format('\u0394r'),  '{:>28.10f}'.format(self.eigvecs[5,i]), '{:>28.12f}'.format(self.phase[5,i]), '{:>28.12f}'.format(self.amps[5,i]))
                    print('{:>28}'.format('\u0394xc'),  '{:>28.10f}'.format(self.eigvecs[6,i]), '{:>28.12f}'.format(self.phase[6,i]), '{:>28.12f}'.format(self.amps[6,i]))
                    print('{:>28}'.format('\u0394yc'),  '{:>28.10f}'.format(self.eigvecs[7,i]), '{:>28.12f}'.format(self.phase[7,i]), '{:>28.12f}'.format(self.amps[7,i]))
                    print('{:>28}'.format('\u0394zc'),  '{:>28.10f}'.format(self.eigvecs[8,i]), '{:>28.12f}'.format(self.phase[8,i]), '{:>28.12f}'.format(self.amps[8,i]))
                    print('{:>28}'.format('\u0394\u03C6:'),  '{:>28.10f}'.format(self.eigvecs[9,i]), '{:>28.12f}'.format(self.phase[9,i]), '{:>28.12f}'.format(self.amps[9,i]))
                    print('{:>28}'.format('\u0394\u03B8:'),  '{:>28.10f}'.format(self.eigvecs[10,i]), '{:>28.12f}'.format(self.phase[10,i]), '{:>28.12f}'.format(self.amps[10,i]))
                    print('{:>28}'.format('\u0394\u03C8:'),  '{:>28.10f}'.format(self.eigvecs[11,i]), '{:>28.12f}'.format(self.phase[11,i]), '{:>28.12f}'.format(self.amps[11,i]))
                else:
                    print('{:>28}'.format('component'),  '{:>28}'.format('eigenvector'), '{:>28}'.format('phase'), '{:>28}'.format('amplitude'))
                    print('{:>28}'.format('\u0394u'),  '{:>28.10f}'.format(self.eigvecs[0,i]), '{:>28.12f}'.format(self.phase[0,i]), '{:>28.12f}'.format(self.amps[0,i]))
                    print('{:>28}'.format('\u0394v'),  '{:>28.10f}'.format(self.eigvecs[1,i]), '{:>28.12f}'.format(self.phase[1,i]), '{:>28.12f}'.format(self.amps[1,i]))
                    print('{:>28}'.format('\u0394w'),  '{:>28.10f}'.format(self.eigvecs[2,i]), '{:>28.12f}'.format(self.phase[2,i]), '{:>28.12f}'.format(self.amps[2,i]))
                    print('{:>28}'.format('\u0394p'),  '{:>28.10f}'.format(self.eigvecs[3,i]), '{:>28.12f}'.format(self.phase[3,i]), '{:>28.12f}'.format(self.amps[3,i]))
                    print('{:>28}'.format('\u0394q'),  '{:>28.10f}'.format(self.eigvecs[4,i]), '{:>28.12f}'.format(self.phase[4,i]), '{:>28.12f}'.format(self.amps[4,i]))
                    print('{:>28}'.format('\u0394r'),  '{:>28.10f}'.format(self.eigvecs[5,i]), '{:>28.12f}'.format(self.phase[5,i]), '{:>28.12f}'.format(self.amps[5,i]))
                    print('{:>28}'.format('\u0394\u03C6:'),  '{:>28.10f}'.format(self.eigvecs[6,i]), '{:>28.12f}'.format(self.phase[6,i]), '{:>28.12f}'.format(self.amps[6,i]))
                    print('{:>28}'.format('\u0394\u03B8:'),  '{:>28.10f}'.format(self.eigvecs[7,i]), '{:>28.12f}'.format(self.phase[7,i]), '{:>28.12f}'.format(self.amps[7,i]))
                    print('{:>28}'.format('\u0394\u03C8:'),  '{:>28.10f}'.format(self.eigvecs[8,i]), '{:>28.12f}'.format(self.phase[8,i]), '{:>28.12f}'.format(self.amps[8,i]))   
        
    def solve_nondim_dynamics_system_hm(self, print_results = True, remove_xyz = False, norm_types = False):
        
        '''
        Solves the LTI system and report eigensolution. Uses the nondimensional
        linearized equations of motion base on the Moulton and Hunsaker nondimensionalization.
        
        Parameters
        -----------
        
        print_results: boolean
            flag for printing eigenvalue and eigenvector solutions
            
        remove_xyz: boolean
            flag for removing the xyz components from the LTI
            
        norm_types: boolean
            flag for normalizing the eigenvector components relative to the seperate
            state variable types (velocities, rotation rates, positions, orientations)

        '''
                                
        V_o = self.V
        u_o = self.eq_velo[0]
        v_o = self.eq_velo[1]
        w_o = self.eq_velo[2]
        
        # nondimensionalize velocities
        mu_o = u_o/V_o
        beta_o = v_o/V_o
        alpha_o = w_o/V_o   
        
        phi_o = (self.eq_euler[0])
        theta_o = (self.eq_euler[1])
             
        '''INTERNAL CALCULATIONS'''
        
        ST_o = np.sin(theta_o)
        CT_o = np.cos(theta_o)
        TT_o = np.tan(theta_o)
        
        SP_o = np.sin(phi_o)
        CP_o = np.cos(phi_o)

        if self.shss!= True:
            '''If not SHSS, used SCT equations to get pqr, these should match the trim solutions.
            Using this form allows the coordinate system assumption to influence these values.'''
            p_o = (-self.g*SP_o*ST_o*CT_o)/(w_o*ST_o + u_o*CP_o*CT_o)
            q_o = (self.g*SP_o*SP_o*CT_o*CT_o)/(w_o*ST_o + u_o*CP_o*CT_o)
            r_o = (self.g*SP_o*CP_o*CT_o*CT_o)/(w_o*ST_o + u_o*CP_o*CT_o)
        else:
            '''If SHSS, set equal to trim solution values. These should always be zero.
            They are not zero if you used the equations above in a SHSS.'''
            p_o = self.eq_rot[0]
            q_o = self.eq_rot[1]
            r_o = self.eq_rot[2]
        
        # nondimensionalize rotation rates
        p_o = p_o*self.V/self.g
        q_o = q_o*self.V/self.g
        r_o = r_o*self.V/self.g
        
        print('\n')
        
        # nondimensionalize inertias
        ell_xzb = self.Ixzb/self.Ixxb
        ell_zxb = self.Ixzb/self.Izzb
        ell_xyb = self.Ixyb/self.Ixxb
        ell_yxb = self.Ixyb/self.Iyyb
        ell_yzb = self.Iyzb/self.Iyyb
        ell_zyb = self.Iyzb/self.Izzb
        
        # nondimensionalize force derivatives
        Kxb_mu = self.Fxb_u*V_o/self.W
        Kxb_beta = self.Fxb_v*V_o/self.W
        Kxb_alpha = self.Fxb_w*V_o/self.W
        
        Kxb_p = self.Fxb_p*self.g/V_o/self.W
        Kxb_q = self.Fxb_q*self.g/V_o/self.W
        Kxb_r = self.Fxb_r*self.g/V_o/self.W
        
        Kyb_mu = self.Fyb_u*V_o/self.W
        Kyb_beta = self.Fyb_v*V_o/self.W
        Kyb_alpha = self.Fyb_w*V_o/self.W
        
        Kyb_p = self.Fyb_p*self.g/V_o/self.W
        Kyb_q = self.Fyb_q*self.g/V_o/self.W
        Kyb_r = self.Fyb_r*self.g/V_o/self.W
    
        Kzb_mu = self.Fzb_u*V_o/self.W
        Kzb_beta = self.Fzb_v*V_o/self.W
        Kzb_alpha = self.Fzb_w*V_o/self.W
        
        Kzb_p = self.Fzb_p*self.g/V_o/self.W
        Kzb_q = self.Fzb_q*self.g/V_o/self.W
        Kzb_r = self.Fzb_r*self.g/V_o/self.W
        
        # nondimensionalize moment derivatives
        Kellb_mu = self.Mxb_u*(V_o**3)/(self.g**2)/self.Ixxb
        Kellb_beta = self.Mxb_v*(V_o**3)/(self.g**2)/self.Ixxb
        Kellb_alpha = self.Mxb_w*(V_o**3)/(self.g**2)/self.Ixxb
        
        Kellb_p = self.Mxb_p*(V_o)/(self.g)/self.Ixxb
        Kellb_q = self.Mxb_q*(V_o)/(self.g)/self.Ixxb
        Kellb_r = self.Mxb_r*(V_o)/(self.g)/self.Ixxb
        
        Kmb_mu = self.Myb_u*(V_o**3)/(self.g**2)/self.Iyyb
        Kmb_beta = self.Myb_v*(V_o**3)/(self.g**2)/self.Iyyb
        Kmb_alpha = self.Myb_w*(V_o**3)/(self.g**2)/self.Iyyb
        
        Kmb_p = self.Myb_p*(V_o)/(self.g)/self.Iyyb
        Kmb_q = self.Myb_q*(V_o)/(self.g)/self.Iyyb
        Kmb_r = self.Myb_r*(V_o)/(self.g)/self.Iyyb
        
        Knb_mu = self.Mzb_u*(V_o**3)/(self.g**2)/self.Izzb
        Knb_beta = self.Mzb_v*(V_o**3)/(self.g**2)/self.Izzb
        Knb_alpha = self.Mzb_w*(V_o**3)/(self.g**2)/self.Izzb
        
        Knb_p = self.Mzb_p*(V_o)/(self.g)/self.Izzb
        Knb_q = self.Mzb_q*(V_o)/(self.g)/self.Izzb
        Knb_r = self.Mzb_r*(V_o)/(self.g)/self.Izzb
        
        Kzb_alpha_hat = self.Fzb_wdot*self.g/self.W
        Kmb_alpha_hat = self.Myb_wdot*(V_o**2)/(self.g)/self.Ixxb
                
        B_matrix = np.array([[                  1, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0,   1,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,   1 - Kzb_alpha_hat,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,          1,   -ell_xyb,   -ell_xzb, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,      -Kmb_alpha_hat,   -ell_yxb,          1,   -ell_yzb, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,   -ell_zxb,   -ell_zyb,          1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
        
        AMxp = Kellb_p + self.Ixzb/self.Ixxb*q_o - self.Ixyb/self.Ixxb*r_o
        AMxq = Kellb_q - self.hzb/self.Ixxb*V_o/self.g + (self.Iyyb - self.Izzb)/self.Ixxb*r_o + 2*self.Iyzb/self.Ixxb*q_o + self.Ixzb/self.Ixxb*p_o
        AMxr = Kellb_r + self.hyb/self.Ixxb*V_o/self.g + (self.Iyyb - self.Izzb)/self.Ixxb*q_o - 2*self.Iyzb/self.Ixxb*r_o - self.Ixyb/self.Ixxb*p_o
        AMyp = Kmb_p + self.hzb/self.Iyyb*V_o/self.g + (self.Izzb - self.Ixxb)/self.Iyyb*r_o - 2*self.Ixzb/self.Iyyb*p_o - self.Iyzb/self.Iyyb*q_o
        AMyq = Kmb_q + self.Ixyb/self.Iyyb*r_o - self.Iyzb/self.Iyyb*p_o
        AMyr = Kmb_r - self.hxb/self.Iyyb*V_o/self.g + (self.Izzb - self.Ixxb)/self.Iyyb*p_o + 2*self.Ixzb/self.Iyyb*r_o + self.Ixyb/self.Iyyb*q_o
        AMzp = Knb_p - self.hyb/self.Izzb*V_o/self.g + (self.Ixxb - self.Iyyb)/self.Izzb*q_o + 2*self.Ixyb/self.Izzb*p_o + self.Iyzb/self.Izzb*r_o
        AMzq = Knb_q + self.hxb/self.Izzb*V_o/self.g + (self.Ixxb - self.Iyyb)/self.Izzb*p_o - 2*self.Ixyb/self.Izzb*q_o - self.Ixzb/self.Izzb*r_o
        AMzr = Knb_r + self.Iyzb/self.Izzb*p_o - self.Ixzb/self.Izzb*q_o
        
        AxP = beta_o*(CP_o*ST_o) - alpha_o*(SP_o*ST_o)
        AxT = -mu_o*ST_o + beta_o*SP_o*CT_o + alpha_o*CP_o*CT_o
        AxO = -beta_o*(CP_o) + alpha_o*(SP_o)
        AyP = -beta_o*(SP_o) - alpha_o*(CP_o)
        AyT = 0.0
        AyO = mu_o*CT_o + beta_o*(SP_o*ST_o) + alpha_o*(CP_o*ST_o)
        AzP = beta_o*CP_o*CT_o - alpha_o*SP_o*CT_o
        AzT = -mu_o*CT_o - beta_o*SP_o*ST_o - alpha_o*CP_o*ST_o
        AzO = 0.0
        
        A441 = 0
        A442 = (q_o*SP_o + r_o*CP_o)*(1/(CT_o*CT_o))
        A443 = (-q_o*SP_o - r_o*CP_o)
        A444 = 0
        A445 = (q_o*SP_o + r_o*CP_o)*(TT_o/CT_o)
        
        A_matrix = np.array([[              Kxb_mu,             Kxb_beta + r_o,            Kxb_alpha - q_o,                Kxb_p,      Kxb_q - alpha_o,       Kxb_r + beta_o,       0.0, 0.0, 0.0,               0.0,             -CT_o, 0.0],
                             [        Kyb_mu - r_o,                   Kyb_beta,            Kyb_alpha + p_o,      Kyb_p + alpha_o,                Kyb_q,         Kyb_r - mu_o,       0.0, 0.0, 0.0,         CP_o*CT_o,        -SP_o*ST_o, 0.0],
                             [        Kzb_mu + q_o,             Kzb_beta - p_o,                  Kzb_alpha,       Kzb_p - beta_o,         Kzb_q + mu_o,                Kzb_r,       0.0, 0.0, 0.0,        -SP_o*CT_o,        -CP_o*ST_o, 0.0],
                             [            Kellb_mu,                 Kellb_beta,                Kellb_alpha,                 AMxp,                 AMxq,                 AMxr,       0.0, 0.0, 0.0,               0.0,               0.0, 0.0],
                             [              Kmb_mu,                   Kmb_beta,                  Kmb_alpha,                 AMyp,                 AMyq,                 AMyr,       0.0, 0.0, 0.0,               0.0,               0.0, 0.0],
                             [              Knb_mu,                   Knb_beta,                  Knb_alpha,                 AMzp,                 AMzq,                 AMzr,       0.0, 0.0, 0.0,               0.0,               0.0, 0.0],
                             [                CT_o,                  SP_o*ST_o,                  CP_o*ST_o,                  0.0,                  0.0,                  0.0,       0.0, 0.0, 0.0,               AxP,               AxT, AxO],
                             [                   0,                       CP_o,                      -SP_o,                  0.0,                  0.0,                  0.0,       0.0, 0.0, 0.0,               AyP,               AyT, AyO],
                             [               -ST_o,                  SP_o*CT_o,                  CP_o*CT_o,                  0.0,                  0.0,                  0.0,       0.0, 0.0, 0.0,               AzP,               AzT, AzO],
                             [                 0.0,                        0.0,                        0.0,                  1.0,            SP_o*TT_o,            CP_o*TT_o,       0.0, 0.0, 0.0,              A441,              A442, 0.0],
                             [                 0.0,                        0.0,                        0.0,                  0.0,                 CP_o,                -SP_o,       0.0, 0.0, 0.0,              A443,               0.0, 0.0],
                             [                 0.0,                        0.0,                        0.0,                  0.0,            SP_o/CT_o,            CP_o/CT_o,       0.0, 0.0, 0.0,              A444,              A445, 0.0]])
       
        # remove xyz components from A and B matrices
        if remove_xyz == True:
            A_matrix = np.delete(A_matrix, [6,7,8], axis=0)
            A_matrix = np.delete(A_matrix, [6,7,8], axis=1)
            
            B_matrix = np.delete(B_matrix, [6,7,8], axis=0)
            B_matrix = np.delete(B_matrix, [6,7,8], axis=1)

        # calculate C matrix
        C_matrix = np.matmul(np.linalg.inv(B_matrix),A_matrix) # turns this into the special eigenproblem?
        self.eigvals, self.eigvecs = eig(C_matrix)
        
        # redimensionalize eigenvalues
        self.eigvals = self.eigvals*self.g/V_o
        
        # print A, B, and C matrices
        print('A-Matrix')
        print('\n'.join([''.join(['{:>16.6f}'.format(item) for item in row]) 
                for row in A_matrix]))
        
        print('\n')
        print('B-Matrix')
        print('\n'.join([''.join(['{:>16.6f}'.format(item) for item in row]) 
                for row in B_matrix]))
        
        print('\n')
        print('C-Matrix')
        print('\n'.join([''.join(['{:>16.6f}'.format(item) for item in row]) 
                for row in C_matrix]))    

        #normalize eigenvectors relative to the largest in each array
        for i in range(len(self.eigvecs[0,:])):
            
            index_max = np.argmax(np.abs(self.eigvecs[:,i]))
            
            '''MULTIPLYING BY CC JUST SHIFTS THE PHASE ANGLES TO BE RELATIVE TO THAT COMPONENT'''
            # index_max = 8 
            
            cc = np.conj(self.eigvecs[index_max,i])
            
            new_vec = cc*self.eigvecs[:,i]
            
            new_vec = new_vec / np.sqrt(np.sum(np.square(np.abs(new_vec))))
            
            self.eigvecs[:,i] = new_vec

        # normalize the eigenvectors within their specific state variable types
        if norm_types == True:
            
            if remove_xyz == True:
                
                for i in range(len(self.eigvecs[:,0])):
                    # velocities
                    self.eigvecs[0:3,i] = self.eigvecs[0:3,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[0:3,i]))))
                    # rotatation rates
                    self.eigvecs[3:6,i] = self.eigvecs[3:6,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[3:6,i]))))
                    # orientation
                    self.eigvecs[6:9,i] = self.eigvecs[6:9,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[6:9,i]))))
            else:
                
                for i in range(len(self.eigvecs[:,0])):
                    # velocities
                    self.eigvecs[0:3,i] = self.eigvecs[0:3,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[0:3,i]))))
                    # rotation rates
                    self.eigvecs[3:6,i] = self.eigvecs[3:6,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[3:6,i]))))
                    # earth fixed positions
                    self.eigvecs[6:9,i] = self.eigvecs[6:9,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[6:9,i]))))
                    # orientation
                    self.eigvecs[9:12,i] = self.eigvecs[9:12,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[9:12,i]))))            

        # sort eigen vectors according to type first (real and imaginary) and then magnitude
        self.eigvals, i_sort = self.sort_eigenvalues_and_indices(self.eigvals)
        # match eigenvectors to sorted eigenvalues
        self.eigvecs = self.eigvecs[:,i_sort]
        # calculate eigenvector amplitudes and phase
        self.amps = np.abs(self.eigvecs)
        self.phase = np.rad2deg(np.arctan2(np.imag(self.eigvecs),np.real(self.eigvecs)))
        
        self.eigreal = np.real(self.eigvals[:])
        self.eigimag = np.imag(self.eigvals[:])
        self.sigma = -np.real(self.eigvals[:])
        self.omegad = np.abs(np.imag(self.eigvals))
        self.period = 2.0*np.pi/self.omegad
                   
        if print_results == True:
            # print eigenvalues
            print('\nEigenvalues') 
            print('\n'.join('{:>32.12f}'.format(item) for item in self.eigvals))
                    
            if self.write_output != False:
                with open(self.output_filename, 'w') as export_handle:
                    for i in range(len(self.eigvals)):
                        print('{:<16.12f}{:<16.12f}'.format(np.real(self.eigvals[i]), np.imag(self.eigvals[i])), file=export_handle)       
            
            # print eigenvalues, related oscillatory properties, and eigenvectors
            for i in range(len(self.eigvecs[0,:])):
                
                print('\n')
                print('{:>24}'.format('Eigenvalue'), '{:<18.12f}'.format(self.eigvals[i]))
                print('{:>24}'.format('Real'), '{:<18.12f}'.format(self.eigreal[i]))
                print('{:>24}'.format('Imaginary'), '{:<18.12f}'.format(self.eigimag[i]))
                print('{:>24}'.format('Period'), '{:<18.12f}'.format(self.period[i]))
                print('{:>24}'.format('Damping'), '{:<18.12f}'.format(self.sigma[i]))
                
                if len(self.eigvecs[:,0]) == 12:
                    print('{:>28}'.format('component'),  '{:>28}'.format('eigenvector'), '{:>28}'.format('phase'), '{:>28}'.format('amplitude'))
                    print('{:>28}'.format('\u0394\u03BC:'), '{:>28.10f}'.format(self.eigvecs[0,i]), '{:>28.12f}'.format(self.phase[0,i]), '{:>28.12f}'.format(self.amps[0,i])) #mu
                    print('{:>28}'.format('\u0394\u03B2:'), '{:>28.10f}'.format(self.eigvecs[1,i]), '{:>28.12f}'.format(self.phase[1,i]), '{:>28.12f}'.format(self.amps[1,i])) #beta
                    print('{:>28}'.format('\u0394\u03B1:'), '{:>28.10f}'.format(self.eigvecs[2,i]), '{:>28.12f}'.format(self.phase[2,i]), '{:>28.12f}'.format(self.amps[2,i])) #alpha
                    print('{:>28}'.format('\u0394pc:'), '{:>28.10f}'.format(self.eigvecs[3,i]), '{:>28.12f}'.format(self.phase[3,i]), '{:>28.12f}'.format(self.amps[3,i])) #phat
                    print('{:>28}'.format('\u0394qc:'), '{:>28.10f}'.format(self.eigvecs[4,i]), '{:>28.12f}'.format(self.phase[4,i]), '{:>28.12f}'.format(self.amps[4,i])) #qhat
                    print('{:>28}'.format('\u0394rc:'), '{:>28.10f}'.format(self.eigvecs[5,i]), '{:>28.12f}'.format(self.phase[5,i]), '{:>28.12f}'.format(self.amps[5,i])) #rhat
                    print('{:>28}'.format('\u0394xc:'), '{:>28.10f}'.format(self.eigvecs[6,i]), '{:>28.12f}'.format(self.phase[6,i]), '{:>28.12f}'.format(self.amps[6,i])) #sigmax
                    print('{:>28}'.format('\u0394yc:'), '{:>28.10f}'.format(self.eigvecs[7,i]), '{:>28.12f}'.format(self.phase[7,i]), '{:>28.12f}'.format(self.amps[7,i])) #sigmay
                    print('{:>28}'.format('\u0394zc:'), '{:>28.10f}'.format(self.eigvecs[8,i]), '{:>28.12f}'.format(self.phase[8,i]), '{:>28.12f}'.format(self.amps[8,i])) #sigmaz
                    print('{:>28}'.format('\u0394\u03C6:'), '{:>28.10f}'.format(self.eigvecs[9,i]), '{:>28.12f}'.format(self.phase[9,i]), '{:>28.12f}'.format(self.amps[9,i])) #PHI
                    print('{:>28}'.format('\u0394\u03B8:'), '{:>28.10f}'.format(self.eigvecs[10,i]), '{:>28.12f}'.format(self.phase[10,i]), '{:>28.12f}'.format(self.amps[10,i])) #THETA       
                    print('{:>28}'.format('\u0394\u03C8:'), '{:>28.10f}'.format(self.eigvecs[11,i]), '{:>28.12f}'.format(self.phase[11,i]), '{:>28.12f}'.format(self.amps[11,i])) #PSI
                else:
                    print('{:>28}'.format('component'),  '{:>28}'.format('eigenvector'), '{:>28}'.format('phase'), '{:>28}'.format('amplitude'))
                    print('{:>28}'.format('\u0394\u03BC:'), '{:>28.10f}'.format(self.eigvecs[0,i]), '{:>28.12f}'.format(self.phase[0,i]), '{:>28.12f}'.format(self.amps[0,i])) #mu
                    print('{:>28}'.format('\u0394\u03B2:'), '{:>28.10f}'.format(self.eigvecs[1,i]), '{:>28.12f}'.format(self.phase[1,i]), '{:>28.12f}'.format(self.amps[1,i])) #beta
                    print('{:>28}'.format('\u0394\u03B1:'), '{:>28.10f}'.format(self.eigvecs[2,i]), '{:>28.12f}'.format(self.phase[2,i]), '{:>28.12f}'.format(self.amps[2,i])) #alpha
                    print('{:>28}'.format('\u0394pc:'), '{:>28.10f}'.format(self.eigvecs[3,i]), '{:>28.12f}'.format(self.phase[3,i]), '{:>28.12f}'.format(self.amps[3,i])) #phat
                    print('{:>28}'.format('\u0394qc:'), '{:>28.10f}'.format(self.eigvecs[4,i]), '{:>28.12f}'.format(self.phase[4,i]), '{:>28.12f}'.format(self.amps[4,i])) #qhat
                    print('{:>28}'.format('\u0394rc:'), '{:>28.10f}'.format(self.eigvecs[5,i]), '{:>28.12f}'.format(self.phase[5,i]), '{:>28.12f}'.format(self.amps[5,i])) #rhat
                    print('{:>28}'.format('\u0394\u03C6:'), '{:>28.10f}'.format(self.eigvecs[6,i]), '{:>28.12f}'.format(self.phase[6,i]), '{:>28.12f}'.format(self.amps[6,i])) #PHI
                    print('{:>28}'.format('\u0394\u03B8:'), '{:>28.10f}'.format(self.eigvecs[7,i]), '{:>28.12f}'.format(self.phase[7,i]), '{:>28.12f}'.format(self.amps[7,i])) #THETA       
                    print('{:>28}'.format('\u0394\u03C8:'), '{:>28.10f}'.format(self.eigvecs[8,i]), '{:>28.12f}'.format(self.phase[8,i]), '{:>28.12f}'.format(self.amps[8,i])) #PSI


    def solve_nondim_dynamics_system_ph(self, print_results = True, remove_xyz = False, norm_types = False):
        
        '''
        Solves the LTI system and report eigensolution. Uses the nondimensional
        linearized equations of motion base on the Phillips nondimensionalization.
        
        Parameters
        -----------
        
        print_results: boolean
            flag for printing eigenvalue and eigenvector solutions
            
        remove_xyz: boolean
            flag for removing the xyz components from the LTI
            
        norm_types: boolean
            flag for normalizing the eigenvector components relative to the seperate
            state variable types (velocities, rotation rates, positions, orientations)

        '''
                        
        V_o = self.V
        u_o = self.eq_velo[0]
        v_o = self.eq_velo[1]
        w_o = self.eq_velo[2]
        
        # nondimensionalize velocities
        mu_o = u_o/V_o
        beta_o = v_o/V_o
        alpha_o = w_o/V_o   
        
        # set reference length
        l_ref = np.sqrt(self.Sw) # square root of wing area
        Ag = l_ref*self.g/(V_o*V_o)
        
        phi_o = (self.eq_euler[0])
        theta_o = (self.eq_euler[1])
        
        '''INTERNAL CALCULATIONS'''
        
        ST_o = np.sin(theta_o)
        CT_o = np.cos(theta_o)
        TT_o = np.tan(theta_o)
        
        SP_o = np.sin(phi_o)
        CP_o = np.cos(phi_o)
        
        if self.shss!= True:
            '''If not SHSS, used SCT equations to get pqr, these should match the trim solutions.
            Using this form allows the coordinate system assumption to influence these values.'''
            p_o = (-self.g*SP_o*ST_o*CT_o)/(w_o*ST_o + u_o*CP_o*CT_o)
            q_o = (self.g*SP_o*SP_o*CT_o*CT_o)/(w_o*ST_o + u_o*CP_o*CT_o)
            r_o = (self.g*SP_o*CP_o*CT_o*CT_o)/(w_o*ST_o + u_o*CP_o*CT_o)
        else:
            '''If SHSS, set equal to trim solution values. These should always be zero.
            They are not zero if you used the equations above in a SHSS.'''
            p_o = self.eq_rot[0]
            q_o = self.eq_rot[1]
            r_o = self.eq_rot[2]
        
        # nondimensionalize rotation rates
        p_o = p_o*l_ref/V_o
        q_o = q_o*l_ref/V_o
        r_o = r_o*l_ref/V_o
        
        print('\n')
        
        W_g = self.W/self.g
                
        # nondimensionalize inertias
        ell_xzb = self.Ixzb/self.Ixxb
        ell_zxb = self.Ixzb/self.Izzb
        ell_xyb = self.Ixyb/self.Ixxb
        ell_yxb = self.Ixyb/self.Iyyb
        ell_yzb = self.Iyzb/self.Iyyb
        ell_zyb = self.Iyzb/self.Izzb
        
        # nondimensionalize force derivatives
        Axb_mu = self.Fxb_u*l_ref/V_o/W_g
        Axb_beta = self.Fxb_v*l_ref/V_o/W_g
        Axb_alpha = self.Fxb_w*l_ref/V_o/W_g
        
        Axb_p = self.Fxb_p/V_o/W_g
        Axb_q = self.Fxb_q/V_o/W_g
        Axb_r = self.Fxb_r/V_o/W_g
        
        Ayb_mu = self.Fyb_u*l_ref/V_o/W_g
        Ayb_beta = self.Fyb_v*l_ref/V_o/W_g
        Ayb_alpha = self.Fyb_w*l_ref/V_o/W_g
        
        Ayb_p = self.Fyb_p/V_o/W_g
        Ayb_q = self.Fyb_q/V_o/W_g
        Ayb_r = self.Fyb_r/V_o/W_g
    
        Azb_mu = self.Fzb_u*l_ref/V_o/W_g
        Azb_beta = self.Fzb_v*l_ref/V_o/W_g
        Azb_alpha = self.Fzb_w*l_ref/V_o/W_g
        
        Azb_p = self.Fzb_p/V_o/W_g
        Azb_q = self.Fzb_q/V_o/W_g
        Azb_r = self.Fzb_r/V_o/W_g
        
        # nondimensionalize moment derivatives
        Aellb_mu = self.Mxb_u*l_ref*l_ref/self.Ixxb/V_o
        Aellb_beta = self.Mxb_v*l_ref*l_ref/self.Ixxb/V_o
        Aellb_alpha = self.Mxb_w*l_ref*l_ref/self.Ixxb/V_o
        
        Aellb_p = self.Mxb_p*l_ref/self.Ixxb/V_o
        Aellb_q = self.Mxb_q*l_ref/self.Ixxb/V_o
        Aellb_r = self.Mxb_r*l_ref/self.Ixxb/V_o
        
        Amb_mu = self.Myb_u*l_ref*l_ref/self.Iyyb/V_o
        Amb_beta = self.Myb_v*l_ref*l_ref/self.Iyyb/V_o
        Amb_alpha = self.Myb_w*l_ref*l_ref/self.Iyyb/V_o
        
        Amb_p = self.Myb_p*l_ref/self.Iyyb/V_o
        Amb_q = self.Myb_q*l_ref/self.Iyyb/V_o
        Amb_r = self.Myb_r*l_ref/self.Iyyb/V_o
        
        Anb_mu = self.Mzb_u*l_ref*l_ref/self.Izzb/V_o
        Anb_beta = self.Mzb_v*l_ref*l_ref/self.Izzb/V_o
        Anb_alpha = self.Mzb_w*l_ref*l_ref/self.Izzb/V_o
        
        Anb_p = self.Mzb_p*l_ref/self.Izzb/V_o
        Anb_q = self.Mzb_q*l_ref/self.Izzb/V_o
        Anb_r = self.Mzb_r*l_ref/self.Izzb/V_o
        
        Azb_alpha_hat = self.Fzb_wdot/W_g
        Amb_alpha_hat = self.Myb_wdot*l_ref/self.Iyyb
        
        # set B matrix
        
        B_matrix = np.array([[                  1, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0,   1,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,   1 - Azb_alpha_hat,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,          1,   -ell_xyb,   -ell_xzb, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,      -Amb_alpha_hat,   -ell_yxb,          1,   -ell_yzb, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,   -ell_zxb,   -ell_zyb,          1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                             [                0.0, 0.0,                 0.0,        0.0,        0.0,        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
        
        AMxp = Aellb_p + self.Ixzb/self.Ixxb*q_o - self.Ixyb/self.Ixxb*r_o
        AMxq = Aellb_q - self.hzb*l_ref/self.Ixxb/V_o + (self.Iyyb - self.Izzb)/self.Ixxb*r_o + 2*self.Iyzb/self.Ixxb*q_o + self.Ixzb/self.Ixxb*p_o
        AMxr = Aellb_r + self.hyb*l_ref/self.Ixxb/V_o + (self.Iyyb - self.Izzb)/self.Ixxb*q_o - 2*self.Iyzb/self.Ixxb*r_o - self.Ixyb/self.Ixxb*p_o
        AMyp = Amb_p + self.hzb*l_ref/self.Iyyb/V_o + (self.Izzb - self.Ixxb)/self.Iyyb*r_o - 2*self.Ixzb/self.Iyyb*p_o - self.Iyzb/self.Iyyb*q_o
        AMyq = Amb_q + self.Ixyb/self.Iyyb*r_o - self.Iyzb/self.Iyyb*p_o
        AMyr = Amb_r - self.hxb*l_ref/self.Iyyb/V_o + (self.Izzb - self.Ixxb)/self.Iyyb*p_o + 2*self.Ixzb/self.Iyyb*r_o + self.Ixyb/self.Iyyb*q_o
        AMzp = Anb_p - self.hyb*l_ref/self.Izzb/V_o + (self.Ixxb - self.Iyyb)/self.Izzb*q_o + 2*self.Ixyb/self.Izzb*p_o + self.Iyzb/self.Izzb*r_o
        AMzq = Anb_q + self.hxb*l_ref/self.Izzb/V_o + (self.Ixxb - self.Iyyb)/self.Izzb*p_o - 2*self.Ixyb/self.Izzb*q_o - self.Ixzb/self.Izzb*r_o
        AMzr = Anb_r + self.Iyzb/self.Izzb*p_o - self.Ixzb/self.Izzb*q_o
        
        AxP = beta_o*(CP_o*ST_o) - alpha_o*(SP_o*ST_o)
        AxT = -mu_o*ST_o + beta_o*SP_o*CT_o + alpha_o*CP_o*CT_o
        AxO = -beta_o*(CP_o) + alpha_o*(SP_o)
        AyP = -beta_o*(SP_o) - alpha_o*(CP_o)
        AyT = 0.0
        AyO = mu_o*CT_o + beta_o*(SP_o*ST_o) + alpha_o*(CP_o*ST_o)
        AzP = beta_o*CP_o*CT_o - alpha_o*SP_o*CT_o
        AzT = -mu_o*CT_o - beta_o*SP_o*ST_o - alpha_o*CP_o*ST_o
        AzO = 0.0

        A441 = 0
        A442 = (q_o*SP_o + r_o*CP_o)*(1/(CT_o*CT_o))
        A443 = (-q_o*SP_o - r_o*CP_o)
        A444 = 0
        A445 = (q_o*SP_o + r_o*CP_o)*(TT_o/CT_o)
        
        # set A matrix
        A_matrix = np.array([[              Axb_mu,             Axb_beta + r_o,            Axb_alpha - q_o,                Axb_p,      Axb_q - alpha_o,       Axb_r + beta_o,       0.0, 0.0, 0.0,               0.0,          -Ag*CT_o, 0.0],
                             [        Ayb_mu - r_o,                   Ayb_beta,            Ayb_alpha + p_o,      Ayb_p + alpha_o,                Ayb_q,         Ayb_r - mu_o,       0.0, 0.0, 0.0,      Ag*CP_o*CT_o,     -Ag*SP_o*ST_o, 0.0],
                             [        Azb_mu + q_o,             Azb_beta - p_o,                  Azb_alpha,       Azb_p - beta_o,         Azb_q + mu_o,                Azb_r,       0.0, 0.0, 0.0,     -Ag*SP_o*CT_o,     -Ag*CP_o*ST_o, 0.0],
                             [            Aellb_mu,                 Aellb_beta,                Aellb_alpha,                 AMxp,                 AMxq,                 AMxr,       0.0, 0.0, 0.0,               0.0,               0.0, 0.0],
                             [              Amb_mu,                   Amb_beta,                  Amb_alpha,                 AMyp,                 AMyq,                 AMyr,       0.0, 0.0, 0.0,               0.0,               0.0, 0.0],
                             [              Anb_mu,                   Anb_beta,                  Anb_alpha,                 AMzp,                 AMzq,                 AMzr,       0.0, 0.0, 0.0,               0.0,               0.0, 0.0],
                             [                CT_o,                  SP_o*ST_o,                  CP_o*ST_o,                  0.0,                  0.0,                  0.0,       0.0, 0.0, 0.0,               AxP,               AxT, AxO],
                             [                   0,                       CP_o,                      -SP_o,                  0.0,                  0.0,                  0.0,       0.0, 0.0, 0.0,               AyP,               AyT, AyO],
                             [               -ST_o,                  SP_o*CT_o,                  CP_o*CT_o,                  0.0,                  0.0,                  0.0,       0.0, 0.0, 0.0,               AzP,               AzT, AzO],
                             [                 0.0,                        0.0,                        0.0,                  1.0,            SP_o*TT_o,            CP_o*TT_o,       0.0, 0.0, 0.0,              A441,              A442, 0.0],
                             [                 0.0,                        0.0,                        0.0,                  0.0,                 CP_o,                -SP_o,       0.0, 0.0, 0.0,              A443,               0.0, 0.0],
                             [                 0.0,                        0.0,                        0.0,                  0.0,            SP_o/CT_o,            CP_o/CT_o,       0.0, 0.0, 0.0,              A444,              A445, 0.0]])


        # remove xyz components from A and B matrices
        if remove_xyz == True:
            A_matrix = np.delete(A_matrix, [6,7,8], axis=0)
            A_matrix = np.delete(A_matrix, [6,7,8], axis=1)
            
            B_matrix = np.delete(B_matrix, [6,7,8], axis=0)
            B_matrix = np.delete(B_matrix, [6,7,8], axis=1)
                    
        # calculate C matrix
        C_matrix = np.matmul(np.linalg.inv(B_matrix),A_matrix)
        self.eigvals, self.eigvecs = eig(C_matrix)
        self.eigvals = self.eigvals*V_o/l_ref # redimensionalize eigenvalues
        
        # print A, B, and C matrices
        print('A-Matrix')
        print('\n'.join([''.join(['{:>16.6f}'.format(item) for item in row]) 
                for row in A_matrix]))
        
        print('\n')
        print('B-Matrix')
        print('\n'.join([''.join(['{:>16.6f}'.format(item) for item in row]) 
                for row in B_matrix]))
        
        print('\n')
        print('C-Matrix')
        print('\n'.join([''.join(['{:>16.6f}'.format(item) for item in row]) 
                for row in C_matrix]))    

        #normalize eigenvectors relative to the largest in each array
        for i in range(len(self.eigvecs[0,:])):
            
            index_max = np.argmax(np.abs(self.eigvecs[:,i]))
            
            '''MULTIPLYING BY CC JUST SHIFTS THE PHASE ANGLES TO BE RELATIVE TO THAT COMPONENT'''
            # index_max = 8 
            
            cc = np.conj(self.eigvecs[index_max,i])
            
            new_vec = cc*self.eigvecs[:,i]
            
            new_vec = new_vec / np.sqrt(np.sum(np.square(np.abs(new_vec))))
            
            self.eigvecs[:,i] = new_vec

        # normalize the eigenvectors within their specific state variable types
        if norm_types == True:
            
            if remove_xyz == True:
                
                for i in range(len(self.eigvecs[:,0])):
                    # velocities
                    self.eigvecs[0:3,i] = self.eigvecs[0:3,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[0:3,i]))))
                    # rotatation rates
                    self.eigvecs[3:6,i] = self.eigvecs[3:6,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[3:6,i]))))
                    # orientation
                    self.eigvecs[6:9,i] = self.eigvecs[6:9,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[6:9,i]))))
            else:
                
                for i in range(len(self.eigvecs[:,0])):
                    # velocities
                    self.eigvecs[0:3,i] = self.eigvecs[0:3,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[0:3,i]))))
                    # rotation rates
                    self.eigvecs[3:6,i] = self.eigvecs[3:6,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[3:6,i]))))
                    # earth fixed positions
                    self.eigvecs[6:9,i] = self.eigvecs[6:9,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[6:9,i]))))
                    # orientation
                    self.eigvecs[9:12,i] = self.eigvecs[9:12,i] / np.sqrt(np.sum(np.square(np.abs(self.eigvecs[9:12,i]))))            

        # sort eigen vectors according to type first (real and imaginary) and then magnitude
        self.eigvals, i_sort = self.sort_eigenvalues_and_indices(self.eigvals)
        # match eigenvectors to sorted eigenvalues
        self.eigvecs = self.eigvecs[:,i_sort]
        # calculate eigenvector amplitudes and phase
        self.amps = np.abs(self.eigvecs)
        self.phase = np.rad2deg(np.arctan2(np.imag(self.eigvecs),np.real(self.eigvecs)))
        
        self.eigreal = np.real(self.eigvals[:])
        self.eigimag = np.imag(self.eigvals[:])
        self.sigma = -np.real(self.eigvals[:])
        self.omegad = np.abs(np.imag(self.eigvals))
        self.period = 2.0*np.pi/self.omegad
                   
        if print_results == True:
            # print eigenvalues
            print('\nEigenvalues') 
            print('\n'.join('{:>32.12f}'.format(item) for item in self.eigvals))
                    
            if self.write_output != False:
                with open(self.output_filename, 'w') as export_handle:
                    for i in range(len(self.eigvals)):
                        print('{:<16.12f}{:<16.12f}'.format(np.real(self.eigvals[i]), np.imag(self.eigvals[i])), file=export_handle)       
            
            # print eigenvalues, related oscillatory properties, and eigenvectors
            for i in range(len(self.eigvecs[0,:])):
                
                print('\n')
                print('{:>24}'.format('Eigenvalue'), '{:<18.12f}'.format(self.eigvals[i]))
                print('{:>24}'.format('Real'), '{:<18.12f}'.format(self.eigreal[i]))
                print('{:>24}'.format('Imaginary'), '{:<18.12f}'.format(self.eigimag[i]))
                print('{:>24}'.format('Period'), '{:<18.12f}'.format(self.period[i]))
                print('{:>24}'.format('Damping'), '{:<18.12f}'.format(self.sigma[i]))
                
                if len(self.eigvecs[:,0]) == 12:
                    print('{:>28}'.format('component'),  '{:>28}'.format('eigenvector'), '{:>28}'.format('phase'), '{:>28}'.format('amplitude'))
                    print('{:>28}'.format('\u0394\u03BC:'), '{:>28.10f}'.format(self.eigvecs[0,i]), '{:>28.12f}'.format(self.phase[0,i]), '{:>28.12f}'.format(self.amps[0,i])) #mu
                    print('{:>28}'.format('\u0394\u03B2:'), '{:>28.10f}'.format(self.eigvecs[1,i]), '{:>28.12f}'.format(self.phase[1,i]), '{:>28.12f}'.format(self.amps[1,i])) #beta
                    print('{:>28}'.format('\u0394\u03B1:'), '{:>28.10f}'.format(self.eigvecs[2,i]), '{:>28.12f}'.format(self.phase[2,i]), '{:>28.12f}'.format(self.amps[2,i])) #alpha
                    print('{:>28}'.format('\u0394pc:'), '{:>28.10f}'.format(self.eigvecs[3,i]), '{:>28.12f}'.format(self.phase[3,i]), '{:>28.12f}'.format(self.amps[3,i])) #phat
                    print('{:>28}'.format('\u0394qc:'), '{:>28.10f}'.format(self.eigvecs[4,i]), '{:>28.12f}'.format(self.phase[4,i]), '{:>28.12f}'.format(self.amps[4,i])) #qhat
                    print('{:>28}'.format('\u0394rc:'), '{:>28.10f}'.format(self.eigvecs[5,i]), '{:>28.12f}'.format(self.phase[5,i]), '{:>28.12f}'.format(self.amps[5,i])) #rhat
                    print('{:>28}'.format('\u0394xc:'), '{:>28.10f}'.format(self.eigvecs[6,i]), '{:>28.12f}'.format(self.phase[6,i]), '{:>28.12f}'.format(self.amps[6,i])) #sigmax
                    print('{:>28}'.format('\u0394yc:'), '{:>28.10f}'.format(self.eigvecs[7,i]), '{:>28.12f}'.format(self.phase[7,i]), '{:>28.12f}'.format(self.amps[7,i])) #sigmay
                    print('{:>28}'.format('\u0394zc:'), '{:>28.10f}'.format(self.eigvecs[8,i]), '{:>28.12f}'.format(self.phase[8,i]), '{:>28.12f}'.format(self.amps[8,i])) #sigmaz
                    print('{:>28}'.format('\u0394\u03C6:'), '{:>28.10f}'.format(self.eigvecs[9,i]), '{:>28.12f}'.format(self.phase[9,i]), '{:>28.12f}'.format(self.amps[9,i])) #PHI
                    print('{:>28}'.format('\u0394\u03B8:'), '{:>28.10f}'.format(self.eigvecs[10,i]), '{:>28.12f}'.format(self.phase[10,i]), '{:>28.12f}'.format(self.amps[10,i])) #THETA       
                    print('{:>28}'.format('\u0394\u03C8:'), '{:>28.10f}'.format(self.eigvecs[11,i]), '{:>28.12f}'.format(self.phase[11,i]), '{:>28.12f}'.format(self.amps[11,i])) #PSI
                else:
                    print('{:>28}'.format('component'),  '{:>28}'.format('eigenvector'), '{:>28}'.format('phase'), '{:>28}'.format('amplitude'))
                    print('{:>28}'.format('\u0394\u03BC:'), '{:>28.10f}'.format(self.eigvecs[0,i]), '{:>28.12f}'.format(self.phase[0,i]), '{:>28.12f}'.format(self.amps[0,i])) #mu
                    print('{:>28}'.format('\u0394\u03B2:'), '{:>28.10f}'.format(self.eigvecs[1,i]), '{:>28.12f}'.format(self.phase[1,i]), '{:>28.12f}'.format(self.amps[1,i])) #beta
                    print('{:>28}'.format('\u0394\u03B1:'), '{:>28.10f}'.format(self.eigvecs[2,i]), '{:>28.12f}'.format(self.phase[2,i]), '{:>28.12f}'.format(self.amps[2,i])) #alpha
                    print('{:>28}'.format('\u0394pc:'), '{:>28.10f}'.format(self.eigvecs[3,i]), '{:>28.12f}'.format(self.phase[3,i]), '{:>28.12f}'.format(self.amps[3,i])) #phat
                    print('{:>28}'.format('\u0394qc:'), '{:>28.10f}'.format(self.eigvecs[4,i]), '{:>28.12f}'.format(self.phase[4,i]), '{:>28.12f}'.format(self.amps[4,i])) #qhat
                    print('{:>28}'.format('\u0394rc:'), '{:>28.10f}'.format(self.eigvecs[5,i]), '{:>28.12f}'.format(self.phase[5,i]), '{:>28.12f}'.format(self.amps[5,i])) #rhat
                    print('{:>28}'.format('\u0394\u03C6:'), '{:>28.10f}'.format(self.eigvecs[6,i]), '{:>28.12f}'.format(self.phase[6,i]), '{:>28.12f}'.format(self.amps[6,i])) #PHI
                    print('{:>28}'.format('\u0394\u03B8:'), '{:>28.10f}'.format(self.eigvecs[7,i]), '{:>28.12f}'.format(self.phase[7,i]), '{:>28.12f}'.format(self.amps[7,i])) #THETA       
                    print('{:>28}'.format('\u0394\u03C8:'), '{:>28.10f}'.format(self.eigvecs[8,i]), '{:>28.12f}'.format(self.phase[8,i]), '{:>28.12f}'.format(self.amps[8,i])) #PSI
    
    def sort_eigenvalues_and_indices(self, eigenvalues):
        
        '''Sorts the eigen values according to type (imaginary and real) and amplitude'''
        
        def amplitude(z):
            return abs(np.real(z))
    
        # Separate eigenvalues into real and imaginary parts
        real_eigenvalues = [(i, val) for i, val in enumerate(eigenvalues) if np.isreal(val)]
        imaginary_eigenvalues = [(i, val) for i, val in enumerate(eigenvalues) if not np.isreal(val)]
        
        # Sort both lists by amplitude, largest first
        real_eigenvalues.sort(key=lambda x: amplitude(x[1]), reverse=False)
        imaginary_eigenvalues.sort(key=lambda x: amplitude(x[1]), reverse=False)
        
        # Combine sorted real and imaginary eigenvalues
        sorted_eigenvalues = [val for _, val in real_eigenvalues] + [val for _, val in imaginary_eigenvalues]
        
        # Combine indices from real and imaginary parts
        sorted_indices = [i for i, _ in real_eigenvalues] + [i for i, _ in imaginary_eigenvalues]
        
        return sorted_eigenvalues, sorted_indices
                
    def plot_eigvals(self):
        
        '''Simple eigenvalue plot function'''
        
        markers = ['1','o','o','>','<','x','s','s']
        
        plt.figure(0, figsize=(5,4))
        plt.grid(visible=True)
        plt.xlabel('Real')
        plt.ylabel('Imaginary')
        
        for i in range(len(self.eigvals[4:])):
            plt.scatter(np.real(self.eigvals[i+4]),np.imag(self.eigvals[i+4]), marker=markers[i], color='k')
        plt.tight_layout()
        plt.show()        
        
        
if __name__ == "__main__":
    
    
    '''INPUTS'''

    V = 355 #ft/s
    gamma = np.deg2rad(0.0) #rad
    phi = np.deg2rad(0) #rad
    H = 24000. #ft
    cg_shift = [0.0, 0.0, 0.0] #ft
    
    SHSS = False
    COMP = False
    STALL = False
    
    '''RUN CASE'''
    case = dynamicAnalysis(write_output=False, output_filename = 'test.txt',
                            shss=SHSS, compressible=COMP, coords_approx=False, derivs_approx=False,
                            stall=STALL, cg_shift=cg_shift)
    
    case.solve_equilibrium_state(V, H, gamma, phi)
    case.solve_derivatives()
    case.generate_control_matrix()
    
    '''SET LTI System type
    TA - dimensional system from dissertation
    hm - nondimensionalized system using Hunsaker/Moulton nondimensionalization
    ph - nondimensionalized system using Phillips nondimensionalization
    '''

    dim_type_select = ['TA', 'hm', 'ph']
    
    dim_type = 1
    
    if dim_type_select[dim_type] == 'TA':
    
        case.solve_dynamics_system(remove_xyz = False)
        case.plot_eigvals()    
        
        real = case.eigreal
        imag = case.eigimag
        
        eig_vecs = case.eigvecs
        omegads = case.omegad
        sigmas = case.sigma
        
    elif dim_type_select[dim_type] == 'hm':
    
        case.solve_nondim_dynamics_system_hm(remove_xyz = False)
        case.plot_eigvals()    
        
        real = case.eigreal
        imag = case.eigimag
        
        eig_vecs = case.eigvecs
        omegads = case.omegad
        sigmas = case.sigma
        
    elif dim_type_select[dim_type] == 'ph':
    
        case.solve_nondim_dynamics_system_ph(remove_xyz = False)
        case.plot_eigvals()    
        
        real = case.eigreal
        imag = case.eigimag
        
        eig_vecs = case.eigvecs
        omegads = case.omegad
        sigmas = case.sigma
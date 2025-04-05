import numpy as np

class derivs:
    '''This is a stand alone class for storing the necessary derivatives for 
    the linearized dynamic analysis code.'''
    def __init__(self):
        self.Fxb_u = 0.0
        self.Fxb_v = 0.0
        self.Fxb_w = 0.0
        self.Fxb_p = 0.0
        self.Fxb_q = 0.0
        self.Fxb_r = 0.0
        
        self.Fyb_u = 0.0
        self.Fyb_v = 0.0
        self.Fyb_w = 0.0
        self.Fyb_p = 0.0
        self.Fyb_q = 0.0
        self.Fyb_r = 0.0
        
        self.Fzb_u = 0.0
        self.Fzb_v = 0.0
        self.Fzb_w = 0.0
        self.Fzb_p = 0.0
        self.Fzb_q = 0.0
        self.Fzb_r = 0.0
                
        self.Mxb_u = 0.0
        self.Mxb_v = 0.0
        self.Mxb_w = 0.0
        self.Mxb_p = 0.0
        self.Mxb_q = 0.0
        self.Mxb_r = 0.0
        
        self.Myb_u = 0.0
        self.Myb_v = 0.0
        self.Myb_w = 0.0
        self.Myb_p = 0.0
        self.Myb_q = 0.0
        self.Myb_r = 0.0
        
        self.Mzb_u = 0.0
        self.Mzb_v = 0.0
        self.Mzb_w = 0.0
        self.Mzb_p = 0.0
        self.Mzb_q = 0.0
        self.Mzb_r = 0.0
        
        self.Fzb_wdot = 0.0
        self.Myb_wdot = 0.0

        self.Fxb_udot = 0.0
        self.Fxb_vdot = 0.0
        self.Fxb_wdot = 0.0
        self.Fyb_udot = 0.0
        self.Fyb_vdot = 0.0
        self.Fyb_wdot = 0.0
        self.Fzb_udot = 0.0
        self.Fzb_vdot = 0.0
        
        self.Myb_udot = 0.0
                
        self.Fxb_da = 0.0
        self.Fxb_de = 0.0  
        self.Fxb_dr = 0.0
        self.Fxb_tau = 0.0
        
        self.Fyb_da = 0.0
        self.Fyb_de = 0.0
        self.Fyb_dr = 0.0
            
        self.Fzb_da = 0.0
        self.Fzb_de = 0.0
        self.Fzb_dr = 0.0
        
        
        self.Mxb_da = 0.0
        self.Mxb_de = 0.0
        self.Mxb_dr = 0.0
        self.Mxb_tau = 0.0

        self.Myb_da = 0.0
        self.Myb_de = 0.0
        self.Myb_dr = 0.0
        self.Myb_tau = 0.0
        
        self.Mzb_da = 0.0
        self.Mzb_de = 0.0
        self.Mzb_dr = 0.0
        self.Mzb_tau = 0.0
        
        self.Ixxb = 0.0
        self.Iyyb = 0.0
        self.Izzb = 0.0
        self.Ixyb = 0.0
        self.Ixzb = 0.0
        self.Iyzb = 0.0

class solveDerivatives:
    
    def __init__(self, aeroModel, aircraft_properties, cg_shift, trim_solution, 
                 compressible, stall, coords_approx = False,
                 derivs_approx = False):
        '''
        
        Estimates and returns the derivatives needed for the linearized analysis.
        
        Parameters
        -----------
        aeroModel: object
            aerodynamic model that will be queried
        aircraft_properties: object
            aircraft property object
        trim_solution: array like
            trim solution object, see trim_functions.py
        cg_shift: array-lie
               vector defining shift in CG
        
        Flags
        ----------
        compressible: boolean
            unused
        stall: boolean
            unused
        coords_approx: boolean
            apply the coordinate system approximation to derivatives and inertia components
        derivs_approx: boolean
            apply the symmetric derivatve approximation
        '''

        self.derivs_sol = derivs() # initialize the derivatives variable for storing solutions
        
        '''Store the necessary input options'''
        self.compressible = compressible
        self.stall = stall
        self.fuselage_FM = False
        self.cg_shift = cg_shift
        self.derivs_approx = derivs_approx
        self.coords_approx = coords_approx
        
        '''Assign the aeroModel, aircraft properties, and thrust model to 
        member variables'''
        self.aeroModel = aeroModel
        self.aircraft_properties = aircraft_properties
        
        # store operating parameters
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
            
        self.hxb = self.aircraft_properties.hx
        self.hyb = self.aircraft_properties.hy
        self.hzb = self.aircraft_properties.hz
        
        # EXTRACT TRIM SOLUTION VALUES
        tau, alpha, beta, da, de, dr = trim_solution.x
        u, v, w, p, q, r, phi, theta = trim_solution.states
        
        self.alpha = alpha
        self.beta = beta
        
        self.eq_velo = np.array([u,v,w])
        self.eq_rot = np.array([p,q,r])
        self.eq_euler = np.array([phi,theta])
        self.eq_inputs = np.array([tau, alpha, beta, da, de, dr])
        
        self.eq_FM = trim_solution.FM_dim
        self.eq_FM_wind = trim_solution.FM

    def solve_derivs(self):
        '''
        Solves for the aerodynamic derivatives at the specified trim condition
        
        Returns
        -----------
        self.derivs_sol: object
            derivative storage object
        '''
        # step sizes used in finite difference estimates
        du = 2 #ft/s
        dv = 1 #ft/s
        dw = 1 #ft/s
        dp = 0.06; #rad/s
        dq = 0.5 * dp;
        dr = 0.5 * dp;
        dda = 1.0 # deg
        dde = 1.0 # deg
        ddr = 1.0 # deg
        dtau = 0.05
        
        # solve derivatives
        self.solve_numeric_derivatives_simple(du, dv, dw, dp, dq, dr, dda, dde, ddr,dtau)
        
        # apply coordinate system approximation
        if self.coords_approx == True:
            self.set_phillips_approx(coords=True,derivs=False)
        
        self.set_deriv_solution()
    
        return self.derivs_sol

    def acceleration_derivatives(self):
        
        # Acceleration derivatives are assumed zero. These could be updated
        # to utilize analytic approximations or values derived from other means.
        
        self.Fzb_wdot = 0.0
        self.Myb_wdot = 0.0
        
        self.Fxb_udot = 0.0
        self.Fxb_vdot = 0.0
        self.Fxb_wdot = 0.0
        self.Fyb_udot = 0.0
        self.Fyb_vdot = 0.0
        self.Fyb_wdot = 0.0
        self.Fzb_vdot = 0.0
        
        self.Fzb_udot = 0.0
        self.Myb_udot = 0.0
        
    def force_derivs(self, Fm2, Fm1, Fp1, Fp2, delta):
        '''fourth order central difference
        
        Parameters
        -----------
        Fm1, Fm2, Fp1, Fp2: float
            function values at steps from equilibrium value
        delta: float
            finite difference step size

        Returns
        -----------
        df_prime: float
            derivative value at equilibrium point
        '''
        df_prime = (-Fp2 + 8*Fp1 - 8*Fm1 + Fm2)/(12*delta)
        return df_prime

    def solve_numeric_derivatives_simple(self, du, dv, dw, dp, dq, dr, dda, dde, ddr, dtau):
            
        """
        Compute the numeric derivatives of body-fixed forces and moments using a fourth-order central difference method.
        
        This method calculates the partial derivatives of the body-fixed forces and moments with respect to the 
        body-fixed velocities (u, v, w), rotation rates (p, q, r), and control inputs (da, de, dr, tau). 
        It evaluates the aerodynamic forces and moments at various perturbations to obtain the derivatives.
    
        Parameters
        ----------
        du : float
            Perturbation step size for the body-fixed forward velocity (u).
        dv : float
            Perturbation step size for the body-fixed lateral velocity (v).
        dw : float
            Perturbation step size for the body-fixed vertical velocity (w).
        dp : float
            Perturbation step size for the body-fixed roll rate (p).
        dq : float
            Perturbation step size for the body-fixed pitch rate (q).
        dr : float
            Perturbation step size for the body-fixed yaw rate (r).
        dda : float
            Perturbation step size for the aileron control input (da) in radians.
        dde : float
            Perturbation step size for the elevator control input (de) in radians.
        ddr : float
            Perturbation step size for the rudder control input (dr) in radians.
        dtau : float
            Perturbation step size for the thrust control input (tau).
    
        Notes
        -----
        - The function uses the equilibrium values of velocity, rotation rates, and control inputs as the base states 
          to calculate the derivatives.
        - The results include the derivatives of the aerodynamic forces (Fx, Fy, Fz) and moments (Mx, My, Mz) with 
          respect to the velocities, rotation rates, and control inputs.
        - The method employs a fourth-order central difference scheme to estimate the derivatives.
        """
        
        # set equilibrium values
        u_o = self.eq_velo[0]
        v_o = self.eq_velo[1]
        w_o = self.eq_velo[2]
        V_o = self.V
        
        p_o = (self.eq_rot[0])
        q_o = (self.eq_rot[1])
        r_o = (self.eq_rot[2])
        
        pbar_o = p_o*self.bw/(2.*self.V)
        qbar_o = q_o*self.cw/(2.*self.V)
        rbar_o = r_o*self.bw/(2.*self.V)
                
        tau_o = self.eq_inputs[0]
        alpha_o = (self.eq_inputs[1])
        beta_o = (self.eq_inputs[2])
        da_o = (self.eq_inputs[3])
        de_o = (self.eq_inputs[4])
        dr_o = (self.eq_inputs[5])
        
        # trim solution, including the throttle setting
        FM_dim = self.aeroModel.aero_CG_offset_results(alpha = alpha_o, beta = beta_o, pbar = pbar_o, qbar = qbar_o, rbar = rbar_o, da = da_o, de = de_o, dr = dr_o, tau = tau_o,
                                                   V = self.V, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift,thrust_off = False,
                                                   aircraft_props = self.aircraft_properties)
        FX0, FY0, FZ0, Mx0, My0, Mz0 = FM_dim
        
        '''duvw data'''
        u_array = np.array([u_o - 2*du, u_o - du, u_o + du, u_o + 2*du])
        v_array = np.array([v_o - 2*dv, v_o - dv, v_o + dv, v_o + 2*dv])
        w_array = np.array([w_o - 2*dw, w_o - dw, w_o + dw, w_o + 2*dw])
        
        FM_du = np.zeros((4,6))
        FM_dv = np.zeros((4,6))
        FM_dw = np.zeros((4,6))
        
        # u,v,w will alter alpha, beta, and V
        # because V is changed, pbar, qbar, and rbar will change
        
        for i in range(len(u_array)):
            
            '''delta u solutions'''
            u = u_array[i]
            V = np.sqrt(u*u + v_o*v_o + w_o*w_o)
            alpha = np.arctan2(w_o,u)
            beta = np.arcsin(v_o/V)
            
            pbar = pbar_o*V_o/V
            qbar = qbar_o*V_o/V
            rbar = rbar_o*V_o/V
            
            FX, FY, FZ, Mx, My, Mz = self.aeroModel.aero_CG_offset_results(alpha = alpha, beta = beta, pbar = pbar, qbar = qbar, rbar = rbar, da = da_o, de = de_o, dr = dr_o, tau = tau_o,
                                                       V = V, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            
            FM_du[i,:] = np.array([FX, FY, FZ, Mx, My, Mz])
            
            '''delta v solutions'''
            v = v_array[i]
            V = np.sqrt(v*v + u_o*u_o + w_o*w_o)
            beta = np.arcsin(v/V)
            
            pbar = pbar_o*V_o/V
            qbar = qbar_o*V_o/V
            rbar = rbar_o*V_o/V
            
            FX, FY, FZ, Mx, My, Mz = self.aeroModel.aero_CG_offset_results(alpha = alpha_o, beta = beta, pbar = pbar, qbar = qbar, rbar = rbar, da = da_o, de = de_o, dr = dr_o, tau = tau_o,
                                                       V = V, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            
            FM_dv[i,:] = np.array([FX, FY, FZ, Mx, My, Mz])
            
            '''delta w solutions'''
            w = w_array[i]
            V = np.sqrt(w*w + v_o*v_o + u_o*u_o)
            alpha = np.arctan2(w,u_o)
            beta = np.arcsin(v_o/V)
            
            pbar = pbar_o*V_o/V
            qbar = qbar_o*V_o/V
            rbar = rbar_o*V_o/V
            
            FX, FY, FZ, Mx, My, Mz = self.aeroModel.aero_CG_offset_results(alpha = alpha, beta = beta, pbar = pbar, qbar = qbar, rbar = rbar, da = da_o, de = de_o, dr = dr_o, tau = tau_o,
                                                       V = V, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            
            FM_dw[i,:] = np.array([FX, FY, FZ, Mx, My, Mz])
        
        '''derivative of body-fixed forces with respect to u'''
        self.Fxb_u = self.force_derivs(FM_du[0,0],FM_du[1,0],FM_du[2,0],FM_du[3,0], du)
        self.Fyb_u = self.force_derivs(FM_du[0,1],FM_du[1,1],FM_du[2,1],FM_du[3,1], du)
        self.Fzb_u = self.force_derivs(FM_du[0,2],FM_du[1,2],FM_du[2,2],FM_du[3,2], du)
        
        '''derivative of body-fixed moments with respect to u'''
        self.Mxb_u = self.force_derivs(FM_du[0,3],FM_du[1,3],FM_du[2,3],FM_du[3,3], du)
        self.Myb_u = self.force_derivs(FM_du[0,4],FM_du[1,4],FM_du[2,4],FM_du[3,4], du)
        self.Mzb_u = self.force_derivs(FM_du[0,5],FM_du[1,5],FM_du[2,5],FM_du[3,5], du)
        
        '''derivative of body-fixed forces with respect to v'''
        self.Fxb_v = self.force_derivs(FM_dv[0,0],FM_dv[1,0],FM_dv[2,0],FM_dv[3,0], dv)
        self.Fyb_v = self.force_derivs(FM_dv[0,1],FM_dv[1,1],FM_dv[2,1],FM_dv[3,1], dv)
        self.Fzb_v = self.force_derivs(FM_dv[0,2],FM_dv[1,2],FM_dv[2,2],FM_dv[3,2], dv)
        
        '''derivative of body-fixed moments with respect to v'''
        self.Mxb_v = self.force_derivs(FM_dv[0,3],FM_dv[1,3],FM_dv[2,3],FM_dv[3,3], dv)
        self.Myb_v = self.force_derivs(FM_dv[0,4],FM_dv[1,4],FM_dv[2,4],FM_dv[3,4], dv)
        self.Mzb_v = self.force_derivs(FM_dv[0,5],FM_dv[1,5],FM_dv[2,5],FM_dv[3,5], dv)
        
        '''derivative of body-fixed forces with respect to w'''
        self.Fxb_w = self.force_derivs(FM_dw[0,0],FM_dw[1,0],FM_dw[2,0],FM_dw[3,0], dw)
        self.Fyb_w = self.force_derivs(FM_dw[0,1],FM_dw[1,1],FM_dw[2,1],FM_dw[3,1], dw)
        self.Fzb_w = self.force_derivs(FM_dw[0,2],FM_dw[1,2],FM_dw[2,2],FM_dw[3,2], dw)
        
        '''derivative of body-fixed moments with respect to w'''
        self.Mxb_w = self.force_derivs(FM_dw[0,3],FM_dw[1,3],FM_dw[2,3],FM_dw[3,3], dw)
        self.Myb_w = self.force_derivs(FM_dw[0,4],FM_dw[1,4],FM_dw[2,4],FM_dw[3,4], dw)
        self.Mzb_w = self.force_derivs(FM_dw[0,5],FM_dw[1,5],FM_dw[2,5],FM_dw[3,5], dw)
        
        '''delta rotation rate Data'''
        dpbar = dp*self.bw/(2.*self.V)
        dqbar = dq*self.cw/(2.*self.V)
        drbar = dr*self.bw/(2.*self.V)
        
        Pbar_array = np.array([pbar_o - 2*dpbar, pbar_o - dpbar, pbar_o + dpbar, pbar_o + 2*dpbar])
        Qbar_array = np.array([qbar_o - 2*dqbar, qbar_o - dqbar, qbar_o + dqbar, qbar_o + 2*dqbar])
        Rbar_array = np.array([rbar_o - 2*drbar, rbar_o - drbar, rbar_o + drbar, rbar_o + 2*drbar])
        
        FM_dPbar = np.zeros((4,6))
        FM_dQbar = np.zeros((4,6))
        FM_dRbar = np.zeros((4,6))
        
        for i in range(len(Pbar_array)):
            pbar = Pbar_array[i]
            qbar = Qbar_array[i]
            rbar = Rbar_array[i]

            
            FXp, FYp, FZp, Mxp, Myp, Mzp = self.aeroModel.aero_CG_offset_results(alpha = alpha_o, beta = beta_o, pbar = pbar, qbar = qbar_o, rbar = rbar_o, da = da_o, de = de_o, dr = dr_o, tau = tau_o,
                                                       V = V_o, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            FM_dPbar[i,:] = np.array([FXp, FYp, FZp, Mxp, Myp, Mzp])
            
            FXq, FYq, FZq, Mxq, Myq, Mzq = self.aeroModel.aero_CG_offset_results(alpha = alpha_o, beta = beta_o, pbar = pbar_o, qbar = qbar, rbar = rbar_o, da = da_o, de = de_o, dr = dr_o, tau = tau_o,
                                                       V = V_o, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            FM_dQbar[i,:] = np.array([FXq, FYq, FZq, Mxq, Myq, Mzq])
            
            FXr, FYr, FZr, Mxr, Myr, Mzr = self.aeroModel.aero_CG_offset_results(alpha = alpha_o, beta = beta_o, pbar = pbar_o, qbar = qbar_o, rbar = rbar, da = da_o, de = de_o, dr = dr_o, tau = tau_o,
                                                       V = V_o, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            FM_dRbar[i,:] = np.array([FXr, FYr, FZr, Mxr, Myr, Mzr])
        
        self.Fxb_p = self.force_derivs(FM_dPbar[0,0],FM_dPbar[1,0],FM_dPbar[2,0],FM_dPbar[3,0], dp)
        self.Fyb_p = self.force_derivs(FM_dPbar[0,1],FM_dPbar[1,1],FM_dPbar[2,1],FM_dPbar[3,1], dp)
        self.Fzb_p = self.force_derivs(FM_dPbar[0,2],FM_dPbar[1,2],FM_dPbar[2,2],FM_dPbar[3,2], dp)
            
        self.Mxb_p = self.force_derivs(FM_dPbar[0,3],FM_dPbar[1,3],FM_dPbar[2,3],FM_dPbar[3,3], dp)
        self.Myb_p = self.force_derivs(FM_dPbar[0,4],FM_dPbar[1,4],FM_dPbar[2,4],FM_dPbar[3,4], dp)
        self.Mzb_p = self.force_derivs(FM_dPbar[0,5],FM_dPbar[1,5],FM_dPbar[2,5],FM_dPbar[3,5], dp)
        
        self.Fxb_q = self.force_derivs(FM_dQbar[0,0],FM_dQbar[1,0],FM_dQbar[2,0],FM_dQbar[3,0], dq)
        self.Fyb_q = self.force_derivs(FM_dQbar[0,1],FM_dQbar[1,1],FM_dQbar[2,1],FM_dQbar[3,1], dq)
        self.Fzb_q = self.force_derivs(FM_dQbar[0,2],FM_dQbar[1,2],FM_dQbar[2,2],FM_dQbar[3,2], dq)
            
        self.Mxb_q = self.force_derivs(FM_dQbar[0,3],FM_dQbar[1,3],FM_dQbar[2,3],FM_dQbar[3,3], dq)
        self.Myb_q = self.force_derivs(FM_dQbar[0,4],FM_dQbar[1,4],FM_dQbar[2,4],FM_dQbar[3,4], dq)
        self.Mzb_q = self.force_derivs(FM_dQbar[0,5],FM_dQbar[1,5],FM_dQbar[2,5],FM_dQbar[3,5], dq)
        
        self.Fxb_r = self.force_derivs(FM_dRbar[0,0],FM_dRbar[1,0],FM_dRbar[2,0],FM_dRbar[3,0], dr)
        self.Fyb_r = self.force_derivs(FM_dRbar[0,1],FM_dRbar[1,1],FM_dRbar[2,1],FM_dRbar[3,1], dr)
        self.Fzb_r = self.force_derivs(FM_dRbar[0,2],FM_dRbar[1,2],FM_dRbar[2,2],FM_dRbar[3,2], dr)
            
        self.Mxb_r = self.force_derivs(FM_dRbar[0,3],FM_dRbar[1,3],FM_dRbar[2,3],FM_dRbar[3,3], dr)
        self.Myb_r = self.force_derivs(FM_dRbar[0,4],FM_dRbar[1,4],FM_dRbar[2,4],FM_dRbar[3,4], dr)
        self.Mzb_r = self.force_derivs(FM_dRbar[0,5],FM_dRbar[1,5],FM_dRbar[2,5],FM_dRbar[3,5], dr)
        
        '''delta control input data'''
        dda = dda*np.pi/180.
        dde = dde*np.pi/180.
        ddr = ddr*np.pi/180.
        
        dda_array = np.array([da_o - 2*dda, da_o - dda, da_o + dda, da_o + 2*dda])
        dde_array = np.array([de_o - 2*dde, de_o - dde, de_o + dde, de_o + 2*dde])
        ddr_array = np.array([dr_o - 2*ddr, dr_o - ddr, dr_o + ddr, dr_o + 2*ddr])
        dtau_array = np.array([tau_o - 2*dtau, tau_o - dtau, tau_o + dtau, tau_o + 2*dtau])
        
        FM_dda = np.zeros((4,6))
        FM_dde = np.zeros((4,6))
        FM_ddr = np.zeros((4,6))
        FM_dtau = np.zeros((4,6))
        
        for i in range(len(Pbar_array)):
            da = dda_array[i]
            de = dde_array[i]
            dr = ddr_array[i]
            tau = dtau_array[i]
        
            
            FXda, FYda, FZda, Mxda, Myda, Mzda = self.aeroModel.aero_CG_offset_results(alpha = alpha_o, beta = beta_o, pbar = pbar_o, qbar = qbar_o, rbar = rbar_o, da = da, de = de_o, dr = dr_o, tau = tau_o,
                                                       V = V_o, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            FM_dda[i,:] = np.array([FXda, FYda, FZda, Mxda, Myda, Mzda])
            
            FXde, FYde, FZde, Mxde, Myde, Mzde = self.aeroModel.aero_CG_offset_results(alpha = alpha_o, beta = beta_o, pbar = pbar_o, qbar = qbar_o, rbar = rbar_o, da = da_o, de = de, dr = dr_o, tau = tau_o,
                                                       V = V_o, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            FM_dde[i,:] = np.array([FXde, FYde, FZde, Mxde, Myde, Mzde])
            
            FXdr, FYdr, FZdr, Mxdr, Mydr, Mzdr = self.aeroModel.aero_CG_offset_results(alpha = alpha_o, beta = beta_o, pbar = pbar_o, qbar = qbar_o, rbar = rbar_o, da = da_o, de = de_o, dr = dr, tau = tau_o,
                                                       V = V_o, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            FM_ddr[i,:] = np.array([FXdr, FYdr, FZdr, Mxdr, Mydr, Mzdr])
            
            FXtau, FYtau, FZtau, Mxtau, Mytau, Mztau = self.aeroModel.aero_CG_offset_results(alpha = alpha_o, beta = beta_o, pbar = pbar_o, qbar = qbar_o, rbar = rbar_o, da = da_o, de = de_o, dr = dr_o, tau = tau,
                                                       V = V_o, H = self.H, rho_0 = self.rho_0, rho = self.rho, cg_shift = self.cg_shift, thrust_off = False,
                                                       aircraft_props = self.aircraft_properties)
            
            FM_dtau[i,:] = np.array([FXtau, FYtau, FZtau, Mxtau, Mytau, Mztau])
        
        self.Fxb_da = self.force_derivs(FM_dda[0,0],FM_dda[1,0],FM_dda[2,0],FM_dda[3,0], dda)
        self.Fyb_da = self.force_derivs(FM_dda[0,1],FM_dda[1,1],FM_dda[2,1],FM_dda[3,1], dda)
        self.Fzb_da = self.force_derivs(FM_dda[0,2],FM_dda[1,2],FM_dda[2,2],FM_dda[3,2], dda)
            
        self.Mxb_da = self.force_derivs(FM_dda[0,3],FM_dda[1,3],FM_dda[2,3],FM_dda[3,3], dda)
        self.Myb_da = self.force_derivs(FM_dda[0,4],FM_dda[1,4],FM_dda[2,4],FM_dda[3,4], dda)
        self.Mzb_da = self.force_derivs(FM_dda[0,5],FM_dda[1,5],FM_dda[2,5],FM_dda[3,5], dda)
        
        self.Fxb_de = self.force_derivs(FM_dde[0,0],FM_dde[1,0],FM_dde[2,0],FM_dde[3,0], dde)
        self.Fyb_de = self.force_derivs(FM_dde[0,1],FM_dde[1,1],FM_dde[2,1],FM_dde[3,1], dde)
        self.Fzb_de = self.force_derivs(FM_dde[0,2],FM_dde[1,2],FM_dde[2,2],FM_dde[3,2], dde)
            
        self.Mxb_de = self.force_derivs(FM_dde[0,3],FM_dde[1,3],FM_dde[2,3],FM_dde[3,3], dde)
        self.Myb_de = self.force_derivs(FM_dde[0,4],FM_dde[1,4],FM_dde[2,4],FM_dde[3,4], dde)
        self.Mzb_de = self.force_derivs(FM_dde[0,5],FM_dde[1,5],FM_dde[2,5],FM_dde[3,5], dde)
        
        self.Fxb_dr = self.force_derivs(FM_ddr[0,0],FM_ddr[1,0],FM_ddr[2,0],FM_ddr[3,0], ddr)
        self.Fyb_dr = self.force_derivs(FM_ddr[0,1],FM_ddr[1,1],FM_ddr[2,1],FM_ddr[3,1], ddr)
        self.Fzb_dr = self.force_derivs(FM_ddr[0,2],FM_ddr[1,2],FM_ddr[2,2],FM_ddr[3,2], ddr)

        self.Mxb_dr = self.force_derivs(FM_ddr[0,3],FM_ddr[1,3],FM_ddr[2,3],FM_ddr[3,3], ddr)
        self.Myb_dr = self.force_derivs(FM_ddr[0,4],FM_ddr[1,4],FM_ddr[2,4],FM_ddr[3,4], ddr)
        self.Mzb_dr = self.force_derivs(FM_ddr[0,5],FM_ddr[1,5],FM_ddr[2,5],FM_ddr[3,5], ddr)
        
        self.Fxb_tau = self.force_derivs(FM_dtau[0,0],FM_dtau[1,0],FM_dtau[2,0],FM_dtau[3,0], dtau)
        self.Fyb_tau = self.force_derivs(FM_dtau[0,1],FM_dtau[1,1],FM_dtau[2,1],FM_dtau[3,1], dtau)
        self.Fzb_tau = self.force_derivs(FM_dtau[0,2],FM_dtau[1,2],FM_dtau[2,2],FM_dtau[3,2], dtau)

        self.Mxb_tau = self.force_derivs(FM_dtau[0,3],FM_dtau[1,3],FM_dtau[2,3],FM_dtau[3,3], dtau)
        self.Myb_tau = self.force_derivs(FM_dtau[0,4],FM_dtau[1,4],FM_dtau[2,4],FM_dtau[3,4], dtau)
        self.Mzb_tau = self.force_derivs(FM_dtau[0,5],FM_dtau[1,5],FM_dtau[2,5],FM_dtau[3,5], dtau)   

        # set phillips approximation
        if self.derivs_approx == True:
            self.set_phillips_approx(coords=False,derivs=True)

        # set acceleration derivatives
        self.acceleration_derivatives()
        
        # thrust derivaties if needed
        # self.Fxb_tau = 0.0
        # self.Mxb_tau = 0.0
        # self.Myb_tau = 0.0
        # self.Mzb_tau = 0.0
        
        print('\n')
                
    def set_phillips_approx(self, coords = False, derivs = False):

        '''
        Sets either the coordinate system assumption or symmetric aircraft assumption
        to derivatives and inertial values
        '''
        
        if coords == True:
            
            # convert inertias and derivatives to wind coordinates
            self.convert_all_bf2wind()
            
            self.eq_velo = np.array([self.V,0.0,0.0])
        
        if derivs == True:
            
            # set derivatives associated with the symmetric aircraft assumption to zero
            
            self.Fxb_v = 0.0
            self.Fyb_u = 0.0
            self.Fyb_w = 0.0
            self.Fzb_v = 0.0
            
            self.Mxb_u = 0.0
            self.Mxb_w = 0.0
            self.Myb_v = 0.0
            self.Mzb_u = 0.0
            self.Mzb_w = 0.0
            
            self.Fxb_p = 0.0
            self.Fxb_r = 0.0
            self.Fyb_q = 0.0
            self.Fzb_p = 0.0
            self.Fzb_r = 0.0
            
            self.Mxb_q = 0.0
            self.Myb_p = 0.0
            self.Myb_r = 0.0
            self.Mzb_q = 0.0

    def convert_all_bf2wind(self):
        
        '''
        Convert derivatives and inertias to wind coordinates.
        '''
        ca = np.cos(self.alpha)
        sa = np.sin(self.alpha)
        cb = np.cos(self.beta)
        sb = np.sin(self.beta)
        
        bdy2wd = np.array([[ca*cb, sb, sa*cb],
                          [-ca*sb, cb, -sa*sb],
                          [-sa, 0.0, ca]])
        
        wd2bdy = np.array([[ca*cb, -ca*sb, -sa],
                          [sb, cb, 0.0],
                          [sa*cb, -sa*sb, ca]])
        
        
        Ibdy = np.array([[self.Ixxb, -self.Ixyb, -self.Ixzb],[-self.Ixyb,  self.Iyyb, -self.Iyzb],[-self.Ixzb, -self.Iyzb,  self.Izzb]])

        Iwd = np.matmul(bdy2wd,np.matmul(Ibdy,wd2bdy))
        self.Ixxb = Iwd[0,0]
        self.Ixyb = -Iwd[0,1]
        self.Ixzb = -Iwd[0,2]
        self.Ixyb = -Iwd[1,0]
        self.Iyyb = Iwd[1,1]
        self.Iyzb = -Iwd[1,2]
        self.Ixzb = -Iwd[2,0]
        self.Iyzb = -Iwd[2,1]
        self.Izzb = Iwd[2,2]
        
        Fbuvw = np.matmul(bdy2wd,np.matmul(np.array([[self.Fxb_u, self.Fxb_v, self.Fxb_w],[self.Fyb_u, self.Fyb_v, self.Fyb_w],[self.Fzb_u, self.Fzb_v, self.Fzb_w]]),wd2bdy))
        self.Fxb_u = Fbuvw[0,0]
        self.Fxb_v = Fbuvw[0,1]
        self.Fxb_w = Fbuvw[0,2]
        self.Fyb_u = Fbuvw[1,0]
        self.Fyb_v = Fbuvw[1,1]
        self.Fyb_w = Fbuvw[1,2]
        self.Fzb_u = Fbuvw[2,0]
        self.Fzb_v = Fbuvw[2,1]
        self.Fzb_w = Fbuvw[2,2]
        Fbpqr = np.matmul(bdy2wd,np.matmul(np.array([[self.Fxb_p, self.Fxb_q, self.Fxb_r],[self.Fyb_p, self.Fyb_q, self.Fyb_r],[self.Fzb_p, self.Fzb_q, self.Fzb_r]]),wd2bdy))
        self.Fxb_p = Fbpqr[0,0]
        self.Fxb_q = Fbpqr[0,1]
        self.Fxb_r = Fbpqr[0,2]
        self.Fyb_p = Fbpqr[1,0]
        self.Fyb_q = Fbpqr[1,1]
        self.Fyb_r = Fbpqr[1,2]
        self.Fzb_p = Fbpqr[2,0]
        self.Fzb_q = Fbpqr[2,1]
        self.Fzb_r = Fbpqr[2,2]
        Fbdadedr = np.matmul(bdy2wd,np.matmul(np.array([[self.Fxb_da, self.Fxb_de, self.Fxb_dr],[self.Fyb_da, self.Fyb_de, self.Fyb_dr],[self.Fzb_da, self.Fzb_de, self.Fzb_dr]]),wd2bdy))
        self.Fxb_da = Fbdadedr[0,0]
        self.Fxb_de = Fbdadedr[0,1]
        self.Fxb_dr = Fbdadedr[0,2]
        self.Fyb_da = Fbdadedr[1,0]
        self.Fyb_de = Fbdadedr[1,1]
        self.Fyb_dr = Fbdadedr[1,2]
        self.Fzb_da = Fbdadedr[2,0]
        self.Fzb_de = Fbdadedr[2,1]
        self.Fzb_dr = Fbdadedr[2,2]
        Fbuvwdot = np.matmul(bdy2wd,np.matmul(np.array([[self.Fxb_udot, self.Fxb_vdot, self.Fxb_wdot],[self.Fyb_udot, self.Fyb_vdot, self.Fyb_wdot],[self.Fzb_udot, self.Fzb_vdot, self.Fzb_wdot]]),wd2bdy))
        self.Fxb_udot = Fbuvwdot[0,0]
        self.Fxb_vdot = Fbuvwdot[0,1]
        self.Fxb_wdot = Fbuvwdot[0,2]
        self.Fyb_udot = Fbuvwdot[1,0]
        self.Fyb_vdot = Fbuvwdot[1,1]
        self.Fyb_wdot = Fbuvwdot[1,2]
        self.Fzb_udot = Fbuvwdot[2,0]
        self.Fzb_vdot = Fbuvwdot[2,1]
        self.Fzb_wdot = Fbuvwdot[2,2]
        
        Mbuvw = np.matmul(bdy2wd,np.matmul(np.array([[self.Mxb_u, self.Mxb_v, self.Mxb_w],[self.Myb_u, self.Myb_v, self.Myb_w],[self.Mzb_u, self.Mzb_v, self.Mzb_w]]),wd2bdy))
        self.Mxb_u = Mbuvw[0,0]
        self.Mxb_v = Mbuvw[0,1]
        self.Mxb_w = Mbuvw[0,2]
        self.Myb_u = Mbuvw[1,0]
        self.Myb_v = Mbuvw[1,1]
        self.Myb_w = Mbuvw[1,2]
        self.Mzb_u = Mbuvw[2,0]
        self.Mzb_v = Mbuvw[2,1]
        self.Mzb_w = Mbuvw[2,2]
        Mbpqr = np.matmul(bdy2wd,np.matmul(np.array([[self.Mxb_p, self.Mxb_q, self.Mxb_r],[self.Myb_p, self.Myb_q, self.Myb_r],[self.Mzb_p, self.Mzb_q, self.Mzb_r]]),wd2bdy))
        self.Mxb_p = Mbpqr[0,0]
        self.Mxb_q = Mbpqr[0,1]
        self.Mxb_r = Mbpqr[0,2]
        self.Myb_p = Mbpqr[1,0]
        self.Myb_q = Mbpqr[1,1]
        self.Myb_r = Mbpqr[1,2]
        self.Mzb_p = Mbpqr[2,0]
        self.Mzb_q = Mbpqr[2,1]
        self.Mzb_r = Mbpqr[2,2]
        Mbdadedr = np.matmul(bdy2wd,np.matmul(np.array([[self.Mxb_da, self.Mxb_de, self.Mxb_dr],[self.Myb_da, self.Myb_de, self.Myb_dr],[self.Mzb_da, self.Mzb_de, self.Mzb_dr]]),wd2bdy))
        self.Mxb_da = Mbdadedr[0,0]
        self.Mxb_de = Mbdadedr[0,1]
        self.Mxb_dr = Mbdadedr[0,2]
        self.Myb_da = Mbdadedr[1,0]
        self.Myb_de = Mbdadedr[1,1]
        self.Myb_dr = Mbdadedr[1,2]
        self.Mzb_da = Mbdadedr[2,0]
        self.Mzb_de = Mbdadedr[2,1]
        self.Mzb_dr = Mbdadedr[2,2]
        Mbuvwdot = np.matmul(bdy2wd,np.matmul(np.array([[0.0, 0.0, 0.0],[self.Myb_udot, 0.0, self.Myb_wdot],[0.0, 0.0, 0.0]]),wd2bdy))
        self.Myb_udot = Mbuvwdot[1,0]
        self.Myb_wdot = Mbuvwdot[1,2]

    def print_derivs(self, print_control = False, print_latex_tables = False):
        
        ''' Print formatted derivatives'''
        
        derivs_mat =    np.array([[self.Fxb_u,self.Fyb_u,self.Fzb_u,self.Mxb_u,self.Myb_u,self.Mzb_u],
                        [self.Fxb_v,self.Fyb_v,self.Fzb_v,self.Mxb_v,self.Myb_v,self.Mzb_v],
                        [self.Fxb_w,self.Fyb_w,self.Fzb_w,self.Mxb_w,self.Myb_w,self.Mzb_w],
                        [self.Fxb_p,self.Fyb_p,self.Fzb_p,self.Mxb_p,self.Myb_p,self.Mzb_p],
                        [self.Fxb_q,self.Fyb_q,self.Fzb_q,self.Mxb_q,self.Myb_q,self.Mzb_q],
                        [self.Fxb_r,self.Fyb_r,self.Fzb_r,self.Mxb_r,self.Myb_r,self.Mzb_r]])
        
        print('Fxb,u, Fyb,u, Fzb,u, Mxb,u, Myb,u, Mzb,u')
        print('Fxb,v, Fyb,v, Fzb,v, Mxb,v, Myb,v, Mzb,v')
        print('Fxb,w, Fyb,w, Fzb,w, Mxb,w, Myb,w, Mzb,w')
        print('Fxb,p, Fyb,p, Fzb,p, Mxb,p, Myb,p, Mzb,p')
        print('Fxb,q, Fyb,q, Fzb,q, Mxb,q, Myb,q, Mzb,q')
        print('Fxb,r, Fyb,r, Fzb,r, Mxb,r, Myb,r, Mzb,r')
        # dstr = '{:>24.16f}{:>24.16f}{:>24.16f}{:>24.16f}{:>24.16f}{:>24.16f}'
        dstr = '{:>24.4f}{:>24.4f}{:>24.4f}{:>24.4f}{:>24.4f}{:>24.4f}'
        print(dstr.format(self.Fxb_u,self.Fyb_u,self.Fzb_u,self.Mxb_u,self.Myb_u,self.Mzb_u))
        print(dstr.format(self.Fxb_v,self.Fyb_v,self.Fzb_v,self.Mxb_v,self.Myb_v,self.Mzb_v))
        print(dstr.format(self.Fxb_w,self.Fyb_w,self.Fzb_w,self.Mxb_w,self.Myb_w,self.Mzb_w))
        print(dstr.format(self.Fxb_p,self.Fyb_p,self.Fzb_p,self.Mxb_p,self.Myb_p,self.Mzb_p))
        print(dstr.format(self.Fxb_q,self.Fyb_q,self.Fzb_q,self.Mxb_q,self.Myb_q,self.Mzb_q))
        print(dstr.format(self.Fxb_r,self.Fyb_r,self.Fzb_r,self.Mxb_r,self.Myb_r,self.Mzb_r))
        print('\n')
        
        if print_control == True:

            print('Fxb,da, Fyb,da, Fzb,da, Mxb,da, Myb,da, Mzb,da')
            print('Fxb,de, Fyb,de, Fzb,de, Mxb,de, Myb,de, Mzb,de')
            print('Fxb,dr, Fyb,dr, Fzb,dr, Mxb,dr, Myb,dr, Mzb,dr')
            print('Fxb,tau, Fyb,tau, Fzb,tau, Mxb,tau, Myb,tau, Mzb,tau')
            dstr = '{:>24.16f}{:>24.16f}{:>24.16f}{:>24.16f}{:>24.16f}{:>24.16f}'
            print(dstr.format(self.Fxb_da,self.Fyb_da,self.Fzb_da,self.Mxb_da,self.Myb_da,self.Mzb_da))
            print(dstr.format(self.Fxb_de,self.Fyb_de,self.Fzb_de,self.Mxb_de,self.Myb_de,self.Mzb_de))
            print(dstr.format(self.Fxb_dr,self.Fyb_dr,self.Fzb_dr,self.Mxb_dr,self.Myb_dr,self.Mzb_dr))
            print(dstr.format(self.Fxb_tau,self.Fyb_tau,self.Fzb_tau,self.Mxb_tau,self.Myb_tau,self.Mzb_tau))
            print('\n')
            
        if print_latex_tables == True:
            
            row_strings = ['$\\bm{u}$','$\\bm{v}$','$\\bm{w}$','$\\bm{p}$','$\\bm{q}$','$\\bm{r}$']
            
            print('\n')
            for i in range(len(derivs_mat[:,0])):
                row = derivs_mat[i,:]

                print('{:^18}'.format(row_strings[i]) + ' & ' + ' & '.join(['{:^16.4f}'.format(item) for item in row]) + '\\\\')  

    def set_deriv_solution(self):
        
        '''Set force derivatives wrt body-fixed velocities and rotation rates'''
        self.derivs_sol.Fxb_u = self.Fxb_u
        self.derivs_sol.Fxb_v = self.Fxb_v
        self.derivs_sol.Fxb_w = self.Fxb_w
        self.derivs_sol.Fxb_p = self.Fxb_p
        self.derivs_sol.Fxb_q = self.Fxb_q
        self.derivs_sol.Fxb_r = self.Fxb_r
        
        self.derivs_sol.Fyb_u = self.Fyb_u
        self.derivs_sol.Fyb_v = self.Fyb_v
        self.derivs_sol.Fyb_w = self.Fyb_w
        self.derivs_sol.Fyb_p = self.Fyb_p
        self.derivs_sol.Fyb_q = self.Fyb_q
        self.derivs_sol.Fyb_r = self.Fyb_r
        
        self.derivs_sol.Fzb_u = self.Fzb_u
        self.derivs_sol.Fzb_v = self.Fzb_v
        self.derivs_sol.Fzb_w = self.Fzb_w
        self.derivs_sol.Fzb_p = self.Fzb_p
        self.derivs_sol.Fzb_q = self.Fzb_q
        self.derivs_sol.Fzb_r = self.Fzb_r

        '''Set moment derivatives wrt body-fixed velocities and rotation rates'''

        self.derivs_sol.Mxb_u = self.Mxb_u
        self.derivs_sol.Mxb_v = self.Mxb_v
        self.derivs_sol.Mxb_w = self.Mxb_w
        self.derivs_sol.Mxb_p = self.Mxb_p
        self.derivs_sol.Mxb_q = self.Mxb_q
        self.derivs_sol.Mxb_r = self.Mxb_r
        
        self.derivs_sol.Myb_u = self.Myb_u
        self.derivs_sol.Myb_v = self.Myb_v
        self.derivs_sol.Myb_w = self.Myb_w
        self.derivs_sol.Myb_p = self.Myb_p
        self.derivs_sol.Myb_q = self.Myb_q
        self.derivs_sol.Myb_r = self.Myb_r
        
        self.derivs_sol.Mzb_u = self.Mzb_u
        self.derivs_sol.Mzb_v = self.Mzb_v
        self.derivs_sol.Mzb_w = self.Mzb_w
        self.derivs_sol.Mzb_p = self.Mzb_p
        self.derivs_sol.Mzb_q = self.Mzb_q
        self.derivs_sol.Mzb_r = self.Mzb_r
        
        '''Set acceleration derivatives'''
        
        self.derivs_sol.Fzb_wdot = self.Fzb_wdot
        self.derivs_sol.Myb_wdot = self.Myb_wdot
        self.derivs_sol.Fxb_udot = self.Fxb_udot
        self.derivs_sol.Fxb_vdot = self.Fxb_vdot
        self.derivs_sol.Fxb_wdot = self.Fxb_wdot
        self.derivs_sol.Fyb_udot = self.Fyb_udot
        self.derivs_sol.Fyb_vdot = self.Fyb_vdot
        self.derivs_sol.Fyb_wdot = self.Fyb_wdot
        self.derivs_sol.Fzb_vdot = self.Fzb_vdot
        
        self.derivs_sol.Fzb_udot = self.Fzb_udot
        self.derivs_sol.Myb_udot = self.Myb_udot
        
        '''SET CONTROL DERIVATIVES'''
        
        self.derivs_sol.Fxb_da = self.Fxb_da
        self.derivs_sol.Fxb_de  = self.Fxb_de 
        self.derivs_sol.Fxb_dr = self.Fxb_dr
        self.derivs_sol.Fxb_tau = self.Fxb_tau
        
        self.derivs_sol.Fyb_da = self.Fyb_da
        self.derivs_sol.Fyb_de = self.Fyb_de
        self.derivs_sol.Fyb_dr = self.Fyb_dr
        self.derivs_sol.Fyb_tau = self.Fyb_tau
        
        self.derivs_sol.Fzb_da = self.Fzb_da
        self.derivs_sol.Fzb_de = self.Fzb_de
        self.derivs_sol.Fzb_dr = self.Fzb_dr
        self.derivs_sol.Fzb_tau = self.Fzb_tau
        
        self.derivs_sol.Mxb_da = self.Mxb_da
        self.derivs_sol.Mxb_de = self.Mxb_de
        self.derivs_sol.Mxb_dr = self.Mxb_dr
        self.derivs_sol.Mxb_tau = self.Mxb_tau

        self.derivs_sol.Myb_da = self.Myb_da
        self.derivs_sol.Myb_de = self.Myb_de
        self.derivs_sol.Myb_dr = self.Myb_dr
        self.derivs_sol.Myb_tau = self.Myb_tau
        
        self.derivs_sol.Mzb_da = self.Mzb_da
        self.derivs_sol.Mzb_de = self.Mzb_de
        self.derivs_sol.Mzb_dr = self.Mzb_dr
        self.derivs_sol.Mzb_tau = self.Mzb_tau
        
        '''Set rotated inertias'''
        self.derivs_sol.Ixxb = self.Ixxb
        self.derivs_sol.Iyyb = self.Iyyb
        self.derivs_sol.Izzb = self.Izzb
        self.derivs_sol.Ixyb = self.Ixyb
        self.derivs_sol.Ixzb = self.Ixzb
        self.derivs_sol.Iyzb = self.Iyzb
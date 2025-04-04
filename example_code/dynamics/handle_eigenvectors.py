import matplotlib.pyplot as plt
import numpy as np

import matplotlib
from matplotlib.ticker import FormatStrFormatter

# f_size = 10
# plt.rcParams.update({'font.size': f_size})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
# plt.rcParams['text.usetex'] = True
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams['axes.axisbelow'] = True


'''Functions useful in running parameter space stability sweeps'''

class derivative_storage:
    
    def __init__(self, nondim_derivs = False):
        self.nondim_derivs = nondim_derivs
        
        self.Fxb_u = []
        self.Fxb_v = []  
        self.Fxb_w = []
        self.Fxb_p = []
        self.Fxb_q = []
        self.Fxb_r = []
        
        self.Fyb_u = []
        self.Fyb_v = []
        self.Fyb_w = []
        self.Fyb_p = []
        self.Fyb_q = []
        self.Fyb_r = []
        
        self.Fzb_u = []
        self.Fzb_v = []
        self.Fzb_w = []
        self.Fzb_p = []
        self.Fzb_q = []
        self.Fzb_r = []
                
        self.Mxb_u = []
        self.Mxb_v = []
        self.Mxb_w = []
        self.Mxb_p = []
        self.Mxb_q = []
        self.Mxb_r = []
        
        self.Myb_u = []
        self.Myb_v = []
        self.Myb_w = []
        self.Myb_p = []
        self.Myb_q = []
        self.Myb_r = []
        
        self.Mzb_u = []
        self.Mzb_v = []
        self.Mzb_w = []
        self.Mzb_p = []
        self.Mzb_q = []
        self.Mzb_r = []
        
        if nondim_derivs == True:
            self.CD_alpha = []
            self.CS_alpha = []
            self.CL_alpha = []
            self.Cl_alpha = []
            self.Cm_alpha = []
            self.Cn_alpha = []
            
            self.CD_beta = []
            self.CS_beta = []
            self.CL_beta = []
            self.Cl_beta = []
            self.Cm_beta = []
            self.Cn_beta = []
            
            self.CD_pbar = []
            self.CS_pbar = []
            self.CL_pbar = []
            self.Cl_pbar = []
            self.Cm_pbar = []
            self.Cn_pbar = []
            
            self.CD_qbar = []
            self.CS_qbar = []
            self.CL_qbar = []
            self.Cl_qbar = []
            self.Cm_qbar = []
            self.Cn_qbar = []
            
            self.CD_rbar = []
            self.CS_rbar = []
            self.CL_rbar = []
            self.Cl_rbar = []
            self.Cm_rbar = []
            self.Cn_rbar = []

                
    def set_nondim_params(self,case):
        
        self.rho = case.rho
        self.V = case.V
        self.Sw = case.Sw
        self.bw = case.bw
        self.cw = case.cw
        
        
    def append_derivatives(self, case):
        
        self.Fxb_u.append(case.Fxb_u)
        self.Fxb_v.append(case.Fxb_v)
        self.Fxb_w.append(case.Fxb_w)
        self.Fxb_p.append(case.Fxb_p)
        self.Fxb_q.append(case.Fxb_q)
        self.Fxb_r.append(case.Fxb_r)
        
        self.Fyb_u.append(case.Fyb_u)
        self.Fyb_v.append(case.Fyb_v)
        self.Fyb_w.append(case.Fyb_w)
        self.Fyb_p.append(case.Fyb_p)
        self.Fyb_q.append(case.Fyb_q)
        self.Fyb_r.append(case.Fyb_r)
        
        self.Fzb_u.append(case.Fzb_u)
        self.Fzb_v.append(case.Fzb_v)
        self.Fzb_w.append(case.Fzb_w)
        self.Fzb_p.append(case.Fzb_p)
        self.Fzb_q.append(case.Fzb_q)
        self.Fzb_r.append(case.Fzb_r)
                
        self.Mxb_u.append(case.Mxb_u)
        self.Mxb_v.append(case.Mxb_v) 
        self.Mxb_w.append(case.Mxb_w)
        self.Mxb_p.append(case.Mxb_p)
        self.Mxb_q.append(case.Mxb_q)
        self.Mxb_r.append(case.Mxb_r) 
        
        self.Myb_u.append(case.Myb_u)
        self.Myb_v.append(case.Myb_v)
        self.Myb_w.append(case.Myb_w)
        self.Myb_p.append(case.Myb_p)
        self.Myb_q.append(case.Myb_q)
        self.Myb_r.append(case.Myb_r)
        
        self.Mzb_u.append(case.Mzb_u)
        self.Mzb_v.append(case.Mzb_v)
        self.Mzb_w.append(case.Mzb_w)
        self.Mzb_p.append(case.Mzb_p)
        self.Mzb_q.append(case.Mzb_q)
        self.Mzb_r.append(case.Mzb_r)

        if self.nondim_derivs == True:
            self.CD_alpha.append(case.deriv_solution.CD_alpha)
            self.CS_alpha.append(case.deriv_solution.CS_alpha)
            self.CL_alpha.append(case.deriv_solution.CL_alpha)
            self.Cl_alpha.append(case.deriv_solution.Cl_alpha)
            self.Cm_alpha.append(case.deriv_solution.Cm_alpha)
            self.Cn_alpha.append(case.deriv_solution.Cn_alpha)
            
            self.CD_beta.append(case.deriv_solution.CD_beta)
            self.CS_beta.append(case.deriv_solution.CS_beta)
            self.CL_beta.append(case.deriv_solution.CL_beta)
            self.Cl_beta.append(case.deriv_solution.Cl_beta)
            self.Cm_beta.append(case.deriv_solution.Cm_beta)
            self.Cn_beta.append(case.deriv_solution.Cn_beta)
            
            self.CD_pbar.append(case.deriv_solution.CD_pbar)
            self.CS_pbar.append(case.deriv_solution.CS_pbar)
            self.CL_pbar.append(case.deriv_solution.CL_pbar)
            self.Cl_pbar.append(case.deriv_solution.Cl_pbar)
            self.Cm_pbar.append(case.deriv_solution.Cm_pbar)
            self.Cn_pbar.append(case.deriv_solution.Cn_pbar)
            
            self.CD_qbar.append(case.deriv_solution.CD_qbar)
            self.CS_qbar.append(case.deriv_solution.CS_qbar)
            self.CL_qbar.append(case.deriv_solution.CL_qbar)
            self.Cl_qbar.append(case.deriv_solution.Cl_qbar)
            self.Cm_qbar.append(case.deriv_solution.Cm_qbar)
            self.Cn_qbar.append(case.deriv_solution.Cn_qbar)
            
            self.CD_rbar.append(case.deriv_solution.CD_rbar)
            self.CS_rbar.append(case.deriv_solution.CS_rbar)
            self.CL_rbar.append(case.deriv_solution.CL_rbar)
            self.Cl_rbar.append(case.deriv_solution.Cl_rbar)
            self.Cm_rbar.append(case.deriv_solution.Cm_rbar)
            self.Cn_rbar.append(case.deriv_solution.Cn_rbar)
        
    def plot_body_fixed_derivatives(self, x_array, baseline_derivatives = None, derivatives_2 = None, case_label='_none',
                           x_label='Bank Angle, $\phi$ [deg]', directory='./',  save_figs=False, num_tics = 5,
                              xlims = (None,None), ylims = [[(None,None), (None,None), (None,None), (None,None)],
                                                            [(None,None), (None,None), (None,None), (None,None)],
                                                            [(None,None), (None,None), (None,None), (None,None)]]):
        
        handletextpad = 0.2
        '''PLOT Body-Fixed Derivatives'''
        figsize = (5.25,5.5)
        # figsize = (6.5,6.81)
        # figsize = (2.925,3.0645)
        
        # X axis
        fig, ax = plt.subplots(2,2,figsize=figsize,sharex=True)

        ax[0,0].plot(x_array,np.array(self.Fxb_u),color='k',linestyle='-', clip_on=False)
        ax[0,0].plot(x_array,np.array(self.Fxb_v),color='k',linestyle='-.', clip_on=False)
        ax[0,0].plot(x_array,np.array(self.Fxb_w),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[0,0].plot(x_array,np.array(baseline_derivatives.Fxb_u),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,0].plot(x_array,np.array(baseline_derivatives.Fxb_v),color='darkgrey',linestyle='-.', clip_on=False)
            ax[0,0].plot(x_array,np.array(baseline_derivatives.Fxb_w),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[0,0].plot(x_array,np.array(derivatives_2.Fxb_u),color='r',linestyle='-', clip_on=False)
            ax[0,0].plot(x_array,np.array(derivatives_2.Fxb_v),color='r',linestyle='-.', clip_on=False)
            ax[0,0].plot(x_array,np.array(derivatives_2.Fxb_w),color='r',linestyle=':', clip_on=False)
        ax[0,0].set_ylabel('Force Derivatives')
        ax[0,0].set_xlim(xlims)
        ax[0,0].set_ylim(ylims[0][0])
        ax[0,0].legend(['$F_{x_b,u}$', '$F_{x_b,v}$', '$F_{x_b,w}$'], handletextpad=handletextpad)
        ax[0,0].grid(visible=True)      
        
        ax[0,1].plot(x_array,np.array(self.Fxb_p),color='k',linestyle='-', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.Fxb_q),color='k',linestyle='-.', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.Fxb_r),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[0,1].plot(x_array,np.array(baseline_derivatives.Fxb_p),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.Fxb_q),color='darkgrey',linestyle='-.', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.Fxb_r),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[0,1].plot(x_array,np.array(derivatives_2.Fxb_p),color='r',linestyle='-', clip_on=False)
            ax[0,1].plot(x_array,np.array(derivatives_2.Fxb_q),color='r',linestyle='-.', clip_on=False)
            ax[0,1].plot(x_array,np.array(derivatives_2.Fxb_r),color='r',linestyle=':', clip_on=False)
        ax[0,1].set_xlim(xlims)
        ax[0,1].set_ylim(ylims[0][1])
        ax[0,1].legend(['$F_{x_b,p}$', '$F_{x_b,q}$', '$F_{x_b,r}$'], handletextpad=handletextpad)
        ax[0,1].grid(visible=True)
        
        ax[1,0].plot(x_array,np.array(self.Mxb_u),color='k',linestyle='-', clip_on=False)
        ax[1,0].plot(x_array,np.array(self.Mxb_v),color='k',linestyle='-.', clip_on=False)
        ax[1,0].plot(x_array,np.array(self.Mxb_w),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Mxb_u),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Mxb_v),color='darkgrey',linestyle='-.', clip_on=False)
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Mxb_w),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[1,0].plot(x_array,np.array(derivatives_2.Mxb_u),color='r',linestyle='-', clip_on=False)
            ax[1,0].plot(x_array,np.array(derivatives_2.Mxb_v),color='r',linestyle='-.', clip_on=False)
            ax[1,0].plot(x_array,np.array(derivatives_2.Mxb_w),color='r',linestyle=':', clip_on=False)
        ax[1,0].set_xlabel(x_label)
        ax[1,0].set_ylabel('Moment Derivatives')
        ax[1,0].set_xlim(xlims)
        ax[1,0].set_ylim(ylims[0][2])
        ax[1,0].legend(['$M_{x_b,u}$', '$M_{x_b,v}$', '$M_{x_b,w}$'], handletextpad=handletextpad)
        ax[1,0].grid(visible=True)
        ax[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        ax[1,1].plot(x_array,np.array(self.Mxb_p),color='k',linestyle='-', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Mxb_q),color='k',linestyle='-.', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Mxb_r),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Mxb_p),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Mxb_q),color='darkgrey',linestyle='-.', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Mxb_r),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[1,1].plot(x_array,np.array(derivatives_2.Mxb_p),color='r',linestyle='-', clip_on=False)
            ax[1,1].plot(x_array,np.array(derivatives_2.Mxb_q),color='r',linestyle='-.', clip_on=False)
            ax[1,1].plot(x_array,np.array(derivatives_2.Mxb_r),color='r',linestyle=':', clip_on=False)
        ax[1,1].set_xlabel(x_label)
        ax[1,1].set_xlim(xlims)
        ax[1,1].set_ylim(ylims[0][3])
        ax[1,1].legend(['$M_{x_b,p}$', '$M_{x_b,q}$', '$M_{x_b,r}$'], handletextpad=handletextpad)
        ax[1,1].grid(visible=True)
        ax[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
    
        plt.tight_layout()
        plt.show()

        if save_figs == True:
            plt.savefig(directory + 'body_X_derivatives_' + case_label + '.pdf')

        # Y axis
        fig, ax = plt.subplots(2,2,figsize=figsize,sharex=True)

        ax[0,0].plot(x_array,np.array(self.Fyb_u),color='k',linestyle='-', clip_on=False)
        ax[0,0].plot(x_array,np.array(self.Fyb_v),color='k',linestyle='-.', clip_on=False)
        ax[0,0].plot(x_array,np.array(self.Fyb_w),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[0,0].plot(x_array,np.array(baseline_derivatives.Fyb_u),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,0].plot(x_array,np.array(baseline_derivatives.Fyb_v),color='darkgrey',linestyle='-.', clip_on=False)
            ax[0,0].plot(x_array,np.array(baseline_derivatives.Fyb_w),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[0,0].plot(x_array,np.array(derivatives_2.Fyb_u),color='r',linestyle='-', clip_on=False)
            ax[0,0].plot(x_array,np.array(derivatives_2.Fyb_v),color='r',linestyle='-.', clip_on=False)
            ax[0,0].plot(x_array,np.array(derivatives_2.Fyb_w),color='r',linestyle=':', clip_on=False)
        ax[0,0].set_ylabel('Force Derivatives')
        ax[0,0].set_xlim(xlims)
        ax[0,0].set_ylim(ylims[1][0])
        ax[0,0].legend(['$F_{y_b,u}$', '$F_{y_b,v}$', '$F_{y_b,w}$'], handletextpad=handletextpad)
        ax[0,0].grid(visible=True)
        
        ax[0,1].plot(x_array,np.array(self.Fyb_p),color='k',linestyle='-', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.Fyb_q),color='k',linestyle='-.', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.Fyb_r),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[0,1].plot(x_array,np.array(baseline_derivatives.Fyb_p),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.Fyb_q),color='darkgrey',linestyle='-.', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.Fyb_r),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[0,1].plot(x_array,np.array(derivatives_2.Fyb_p),color='r',linestyle='-', clip_on=False)
            ax[0,1].plot(x_array,np.array(derivatives_2.Fyb_q),color='r',linestyle='-.', clip_on=False)
            ax[0,1].plot(x_array,np.array(derivatives_2.Fyb_r),color='r',linestyle=':', clip_on=False)
        ax[0,1].set_xlim(xlims)
        ax[0,1].set_ylim(ylims[1][1])
        ax[0,1].legend(['$F_{y_b,p}$', '$F_{y_b,q}$', '$F_{y_b,r}$'], handletextpad=handletextpad)
        ax[0,1].grid(visible=True)
        
        ax[1,0].plot(x_array,np.array(self.Myb_u),color='k',linestyle='-', clip_on=False)
        ax[1,0].plot(x_array,np.array(self.Myb_v),color='k',linestyle='-.', clip_on=False)
        ax[1,0].plot(x_array,np.array(self.Myb_w),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Myb_u),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Myb_v),color='darkgrey',linestyle='-.', clip_on=False)
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Myb_w),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[1,0].plot(x_array,np.array(derivatives_2.Myb_u),color='r',linestyle='-', clip_on=False)
            ax[1,0].plot(x_array,np.array(derivatives_2.Myb_v),color='r',linestyle='-.', clip_on=False)
            ax[1,0].plot(x_array,np.array(derivatives_2.Myb_w),color='r',linestyle=':', clip_on=False)
        ax[1,0].set_xlabel(x_label)
        ax[1,0].set_ylabel('Moment Derivatives')
        ax[1,0].set_xlim(xlims)
        ax[1,0].set_ylim(ylims[1][2])
        ax[1,0].legend(['$M_{y_b,u}$', '$M_{y_b,v}$', '$M_{y_b,w}$'], handletextpad=handletextpad)
        ax[1,0].grid(visible=True)
        ax[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        ax[1,1].plot(x_array,np.array(self.Myb_p),color='k',linestyle='-', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Myb_q),color='k',linestyle='-.', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Myb_r),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Myb_p),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Myb_q),color='darkgrey',linestyle='-.', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Myb_r),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[1,1].plot(x_array,np.array(derivatives_2.Myb_p),color='r',linestyle='-', clip_on=False)
            ax[1,1].plot(x_array,np.array(derivatives_2.Myb_q),color='r',linestyle='-.', clip_on=False)
            ax[1,1].plot(x_array,np.array(derivatives_2.Myb_r),color='r',linestyle=':', clip_on=False)
        ax[1,1].set_xlabel(x_label)
        ax[1,1].set_xlim(xlims)
        ax[1,1].set_ylim(ylims[1][3])
        ax[1,1].legend(['$M_{y_b,p}$', '$M_{y_b,q}$', '$M_{y_b,r}$'], handletextpad=handletextpad)
        ax[1,1].grid(visible=True)
        ax[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        plt.tight_layout()
        plt.show()

        if save_figs == True:
            plt.savefig(directory + 'body_Y_derivatives_' + case_label + '.pdf')
            
        #Z axis
        fig, ax = plt.subplots(2,2,figsize=figsize,sharex=True)
        ax[0,0].plot(x_array,np.array(self.Fzb_u),color='k',linestyle='-', clip_on=False)
        ax[0,0].plot(x_array,np.array(self.Fzb_v),color='k',linestyle='-.', clip_on=False)
        ax[0,0].plot(x_array,np.array(self.Fzb_w),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[0,0].plot(x_array,np.array(baseline_derivatives.Fzb_u),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,0].plot(x_array,np.array(baseline_derivatives.Fzb_v),color='darkgrey',linestyle='-.', clip_on=False)
            ax[0,0].plot(x_array,np.array(baseline_derivatives.Fzb_w),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[0,0].plot(x_array,np.array(derivatives_2.Fzb_u),color='r',linestyle='-', clip_on=False)
            ax[0,0].plot(x_array,np.array(derivatives_2.Fzb_v),color='r',linestyle='-.', clip_on=False)
            ax[0,0].plot(x_array,np.array(derivatives_2.Fzb_w),color='r',linestyle=':', clip_on=False)
        ax[0,0].set_ylabel('Force Derivatives')
        ax[0,0].set_xlim(xlims)
        ax[0,0].set_ylim(ylims[2][0])
        ax[0,0].legend(['$F_{z_b,u}$', '$F_{z_b,v}$', '$F_{z_b,w}$'], handletextpad=handletextpad)
        ax[0,0].grid(visible=True)
        
        ax[0,1].plot(x_array,np.array(self.Fzb_p),color='k',linestyle='-', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.Fzb_q),color='k',linestyle='-.', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.Fzb_r),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[0,1].plot(x_array,np.array(baseline_derivatives.Fzb_p),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.Fzb_q),color='darkgrey',linestyle='-.', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.Fzb_r),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[0,1].plot(x_array,np.array(derivatives_2.Fzb_p),color='r',linestyle='-', clip_on=False)
            ax[0,1].plot(x_array,np.array(derivatives_2.Fzb_q),color='r',linestyle='-.', clip_on=False)
            ax[0,1].plot(x_array,np.array(derivatives_2.Fzb_r),color='r',linestyle=':', clip_on=False)
        ax[0,1].set_xlim(xlims)
        ax[0,1].set_ylim(ylims[2][1])
        ax[0,1].legend(['$F_{z_b,p}$', '$F_{z_b,q}$', '$F_{z_b,r}$'], handletextpad=handletextpad)
        ax[0,1].grid(visible=True)
        
        ax[1,0].plot(x_array,np.array(self.Mzb_u),color='k',linestyle='-', clip_on=False)
        ax[1,0].plot(x_array,np.array(self.Mzb_v),color='k',linestyle='-.', clip_on=False)
        ax[1,0].plot(x_array,np.array(self.Mzb_w),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Mzb_u),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Mzb_v),color='darkgrey',linestyle='-.', clip_on=False)
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Mzb_w),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[1,0].plot(x_array,np.array(derivatives_2.Mzb_u),color='r',linestyle='-', clip_on=False)
            ax[1,0].plot(x_array,np.array(derivatives_2.Mzb_v),color='r',linestyle='-.', clip_on=False)
            ax[1,0].plot(x_array,np.array(derivatives_2.Mzb_w),color='r',linestyle=':', clip_on=False)
        ax[1,0].set_xlabel(x_label)
        ax[1,0].set_ylabel('Moment Derivatives')
        ax[1,0].set_xlim(xlims)
        ax[1,0].set_ylim(ylims[2][2])
        ax[1,0].legend(['$M_{z_b,u}$', '$M_{z_b,v}$', '$M_{z_b,w}$'], handletextpad=handletextpad)
        ax[1,0].grid(visible=True)
        ax[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        ax[1,1].plot(x_array,np.array(self.Mzb_p),color='k',linestyle='-', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Mzb_q),color='k',linestyle='-.', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Mzb_r),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Mzb_p),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Mzb_q),color='darkgrey',linestyle='-.', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Mzb_r),color='darkgrey',linestyle=':', clip_on=False)
        if derivatives_2 != None:
            ax[1,1].plot(x_array,np.array(derivatives_2.Mzb_p),color='r',linestyle='-', clip_on=False)
            ax[1,1].plot(x_array,np.array(derivatives_2.Mzb_q),color='r',linestyle='-.', clip_on=False)
            ax[1,1].plot(x_array,np.array(derivatives_2.Mzb_r),color='r',linestyle=':', clip_on=False)
        ax[1,1].set_xlabel(x_label)
        ax[1,1].set_xlim(xlims)
        ax[1,1].set_ylim(ylims[2][3])
        ax[1,1].legend(['$M_{z_b,p}$', '$M_{z_b,q}$', '$M_{z_b,r}$'], handletextpad=handletextpad)
        ax[1,1].grid(visible=True)
        ax[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        plt.tight_layout()
        plt.show()
        
        if save_figs == True:
            plt.savefig(directory + 'body_Z_derivatives_' + case_label + '.pdf')
        
        
    def plot_nondim_wind_derivatives(self, x_array, baseline_derivatives = None, case_label='_none',
                           x_label='Bank Angle, $\phi$ [deg]', directory='./',  save_figs=False, num_tics = 5,
                              xlims = (None,None), ylims = [[(None,None), (None,None), (None,None), (None,None)],
                                                            [(None,None), (None,None), (None,None), (None,None)],
                                                            [(None,None), (None,None), (None,None), (None,None)]]):
        
        handletextpad = 1.0
        
        '''PLOT Wind Derivatives'''
        figsize = (5.25,5.5)
        
        # X axis
        fig, ax = plt.subplots(2,2,figsize=figsize,sharex=True)
        ax[0,0].plot(x_array,np.array(self.CD_alpha),color='k',linestyle='-', clip_on=False)
        ax[0,0].plot(x_array,np.array(self.CD_beta),color='k',linestyle='-.', clip_on=False)
        if baseline_derivatives != None:
            ax[0,0].plot(x_array,np.array(baseline_derivatives.CD_alpha),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,0].plot(x_array,np.array(baseline_derivatives.CD_beta),color='darkgrey',linestyle='-.', clip_on=False)
        ax[0,0].set_ylabel('Force Coefficient Derivatives')
        ax[0,0].set_xlim(xlims)
        ax[0,0].set_ylim(ylims[0][0])
        ax[0,0].legend(['$C_{D,\\alpha}$', '$C_{D,\\beta}$'], handletextpad=handletextpad)
        ax[0,0].grid(visible=True)
        
        ax[0,1].plot(x_array,np.array(self.CD_pbar),color='k',linestyle='-', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.CD_qbar),color='k',linestyle='-.', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.CD_rbar),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[0,1].plot(x_array,np.array(baseline_derivatives.CD_pbar),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.CD_qbar),color='darkgrey',linestyle='-.', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.CD_rbar),color='darkgrey',linestyle=':', clip_on=False)
        ax[0,1].set_xlim(xlims)
        ax[0,1].set_ylim(ylims[0][1])
        ax[0,1].legend(['$C_{D,\overline{p}}$', '$C_{D,\overline{q}}$', '$C_{D,\overline{r}}$'], handletextpad=handletextpad)
        ax[0,1].grid(visible=True)
        
        ax[1,0].plot(x_array,np.array(self.Cl_alpha),color='k',linestyle='-', clip_on=False)
        ax[1,0].plot(x_array,np.array(self.Cl_beta),color='k',linestyle='-.', clip_on=False)
        if baseline_derivatives != None:
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Cl_alpha),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Cl_beta),color='darkgrey',linestyle='-.', clip_on=False)
        ax[1,0].set_xlabel(x_label)
        ax[1,0].set_ylabel('Moment Coefficient Derivatives')
        ax[1,0].set_xlim(xlims)
        ax[1,0].set_ylim(ylims[0][2])
        ax[1,0].legend(['$C_{\ell,\\alpha}$', '$C_{\ell,\\beta}$'], handletextpad=handletextpad)
        ax[1,0].grid(visible=True)
        ax[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        ax[1,1].plot(x_array,np.array(self.Cl_pbar),color='k',linestyle='-', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Cl_qbar),color='k',linestyle='-.', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Cl_rbar),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Cl_pbar),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Cl_qbar),color='darkgrey',linestyle='-.', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Cl_rbar),color='darkgrey',linestyle=':', clip_on=False)
        ax[1,1].set_xlabel(x_label)
        ax[1,1].set_xlim(xlims)
        ax[1,1].set_ylim(ylims[0][3])
        ax[1,1].legend(['$C_{\ell,\overline{p}}$', '$C_{\ell,\overline{q}}$', '$C_{\ell,\overline{r}}$'], handletextpad=handletextpad)
        ax[1,1].grid(visible=True)
        ax[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        plt.tight_layout()
        plt.show()
        
        if save_figs == True:
            plt.savefig(directory + 'wind_X_derivatives_' + case_label + '.pdf')
        
        # Y axis
        fig, ax = plt.subplots(2,2,figsize=figsize)
        
        ax[0,0].plot(x_array,np.array(self.CS_alpha),color='k',linestyle='-', clip_on=False)
        ax[0,0].plot(x_array,np.array(self.CS_beta),color='k',linestyle='-.', clip_on=False)
        if baseline_derivatives != None:
            ax[0,0].plot(x_array,np.array(baseline_derivatives.CS_alpha),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,0].plot(x_array,np.array(baseline_derivatives.CS_beta),color='darkgrey',linestyle='-.', clip_on=False)
        ax[0,0].set_ylabel('Force Coefficient Derivatives')
        ax[0,0].set_xlim(xlims)
        ax[0,0].set_ylim(ylims[1][0])
        ax[0,0].legend(['$C_{S,\\alpha}$', '$C_{S,\\beta}$'], handletextpad=handletextpad)
        ax[0,0].grid(visible=True)
        
        ax[0,1].plot(x_array,np.array(self.CS_pbar),color='k',linestyle='-', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.CS_qbar),color='k',linestyle='-.', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.CS_rbar),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[0,1].plot(x_array,np.array(baseline_derivatives.CS_pbar),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.CS_qbar),color='darkgrey',linestyle='-.', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.CS_rbar),color='darkgrey',linestyle=':', clip_on=False)
        ax[0,1].set_xlim(xlims)
        ax[0,1].set_ylim(ylims[1][1])
        ax[0,1].legend(['$C_{S,\overline{p}}$', '$C_{S,\overline{q}}$', '$C_{S,\overline{r}}$'], handletextpad=handletextpad)
        ax[0,1].grid(visible=True)
        
        ax[1,0].plot(x_array,np.array(self.Cm_alpha),color='k',linestyle='-', clip_on=False)
        ax[1,0].plot(x_array,np.array(self.Cm_beta),color='k',linestyle='-.', clip_on=False)
        if baseline_derivatives != None:
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Cm_alpha),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Cm_beta),color='darkgrey',linestyle='-.', clip_on=False)
        ax[1,0].set_xlabel(x_label)
        ax[1,0].set_ylabel('Moment Coefficient Derivatives')
        ax[1,0].set_xlim(xlims)
        ax[1,0].set_ylim(ylims[1][2])
        ax[1,0].legend(['$C_{m,\\alpha}$', '$C_{m,\\beta}$'], handletextpad=handletextpad)
        ax[1,0].grid(visible=True)
        ax[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        ax[1,1].plot(x_array,np.array(self.Cm_pbar),color='k',linestyle='-', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Cm_qbar),color='k',linestyle='-.', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Cm_rbar),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Cm_pbar),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Cm_qbar),color='darkgrey',linestyle='-.', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Cm_rbar),color='darkgrey',linestyle=':', clip_on=False)
        ax[1,1].set_xlabel(x_label)
        ax[1,1].set_xlim(xlims)
        ax[1,1].set_ylim(ylims[1][3])
        ax[1,1].legend(['$C_{m,\overline{p}}$', '$C_{m,\overline{q}}$', '$C_{m,\overline{r}}$'], handletextpad=handletextpad)
        ax[1,1].grid(visible=True)
        ax[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        plt.tight_layout()
        plt.show()
        
        if save_figs == True:
            plt.savefig(directory + 'wind_Y_derivatives_' + case_label + '.pdf')
        
        # Z axis
        fig, ax = plt.subplots(2,2,figsize=figsize)
        
        ax[0,0].plot(x_array,np.array(self.CL_alpha),color='k',linestyle='-', clip_on=False)
        ax[0,0].plot(x_array,np.array(self.CL_beta),color='k',linestyle='-.', clip_on=False)
        if baseline_derivatives != None:
            ax[0,0].plot(x_array,np.array(baseline_derivatives.CL_alpha),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,0].plot(x_array,np.array(baseline_derivatives.CL_beta),color='darkgrey',linestyle='-.', clip_on=False)
        ax[0,0].set_ylabel('Force Coefficient Derivatives')
        ax[0,0].set_xlim(xlims)
        ax[0,0].set_ylim(ylims[2][0])
        ax[0,0].legend(['$C_{L,\\alpha}$', '$C_{L,\\beta}$'], handletextpad=handletextpad)
        ax[0,0].grid(visible=True)
        
        ax[0,1].plot(x_array,np.array(self.CL_pbar),color='k',linestyle='-', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.CL_qbar),color='k',linestyle='-.', clip_on=False)
        ax[0,1].plot(x_array,np.array(self.CL_rbar),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[0,1].plot(x_array,np.array(baseline_derivatives.CL_pbar),color='darkgrey',linestyle='-', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.CL_qbar),color='darkgrey',linestyle='-.', clip_on=False)
            ax[0,1].plot(x_array,np.array(baseline_derivatives.CL_rbar),color='darkgrey',linestyle=':', clip_on=False)
        ax[0,1].set_xlim(xlims)
        ax[0,1].set_ylim(ylims[2][1])
        ax[0,1].legend(['$C_{L,\overline{p}}$', '$C_{L,\overline{q}}$', '$C_{L,\overline{r}}$'], handletextpad=handletextpad)
        ax[0,1].grid(visible=True)
        
        ax[1,0].plot(x_array,np.array(self.Cn_alpha),color='k',linestyle='-', clip_on=False)
        ax[1,0].plot(x_array,np.array(self.Cn_beta),color='k',linestyle='-.', clip_on=False)
        if baseline_derivatives != None:
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Cn_alpha),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,0].plot(x_array,np.array(baseline_derivatives.Cn_beta),color='darkgrey',linestyle='-.', clip_on=False)
        # ax[1,0].set_xlabel(x_label)
        ax[1,0].set_ylabel('Moment Coefficient Derivatives')
        ax[1,0].set_xlim(xlims)
        ax[1,0].set_ylim(ylims[0][2])
        ax[1,0].legend(['$C_{n,\\alpha}$', '$C_{n,\\beta}$'], handletextpad=1.0)
        ax[1,0].grid(visible=True)
        ax[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        ax[1,1].plot(x_array,np.array(self.Cn_pbar),color='k',linestyle='-', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Cn_qbar),color='k',linestyle='-.', clip_on=False)
        ax[1,1].plot(x_array,np.array(self.Cn_rbar),color='k',linestyle=':', clip_on=False)
        if baseline_derivatives != None:
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Cn_pbar),color='darkgrey',linestyle='-', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Cn_qbar),color='darkgrey',linestyle='-.', clip_on=False)
            ax[1,1].plot(x_array,np.array(baseline_derivatives.Cn_rbar),color='darkgrey',linestyle=':', clip_on=False)
        ax[1,1].set_xlabel(x_label)
        ax[1,1].set_xlim(xlims)
        ax[1,1].set_ylim(ylims[0][3])
        ax[1,1].legend(['$C_{n,\overline{p}}$', '$C_{n,\overline{q}}$', '$C_{n,\overline{r}}$'], handletextpad=handletextpad)
        ax[1,1].grid(visible=True)
        ax[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        
        plt.tight_layout()
        plt.show()

        if save_figs == True:
            plt.savefig(directory + 'wind_Z_derivatives_' + case_label + '.pdf')

class parse_eigenvectors:
    
    def __init__(self,amps,phase,mode_index,adjust_phase=False):
        
        self.u_amps = []
        self.v_amps = []
        self.w_amps = []
        
        self.p_amps = []
        self.q_amps = []
        self.r_amps = []
        
        self.x_amps = []
        self.y_amps = []
        self.z_amps = []
        
        self.phi_amps = []
        self.theta_amps = []
        self.psi_amps = []
        
        self.u_phase = []
        self.v_phase = []
        self.w_phase = []
        
        self.p_phase = []
        self.q_phase = []
        self.r_phase = []
        
        self.x_phase = []
        self.y_phase = []
        self.z_phase = []
        
        self.phi_phase = []
        self.theta_phase = []
        self.psi_phase = []
        
        # print('Amps: ', len(amps[:]))
        
        for i in range(len(amps[:])):
            self.u_amps.append(amps[i][0,mode_index])
            self.v_amps.append(amps[i][1,mode_index])
            self.w_amps.append(amps[i][2,mode_index])
            
            self.p_amps.append(amps[i][3,mode_index])
            self.q_amps.append(amps[i][4,mode_index])
            self.r_amps.append(amps[i][5,mode_index])
            
            self.x_amps.append(amps[i][6,mode_index])
            self.y_amps.append(amps[i][7,mode_index])
            self.z_amps.append(amps[i][8,mode_index])
            
            self.phi_amps.append(amps[i][9,mode_index])
            self.theta_amps.append(amps[i][10,mode_index])
            self.psi_amps.append(amps[i][11,mode_index])
            
            if adjust_phase != True:
            
                self.u_phase.append(phase[i][0,mode_index])
                self.v_phase.append(phase[i][1,mode_index])
                self.w_phase.append(phase[i][2,mode_index])
                
                self.p_phase.append(phase[i][3,mode_index])
                self.q_phase.append(phase[i][4,mode_index])
                self.r_phase.append(phase[i][5,mode_index])
                
                self.x_phase.append(phase[i][6,mode_index])
                self.y_phase.append(phase[i][7,mode_index])
                self.z_phase.append(phase[i][8,mode_index])
                
                self.phi_phase.append(phase[i][9,mode_index])
                self.theta_phase.append(phase[i][10,mode_index])
                self.psi_phase.append(phase[i][11,mode_index])
            else:
                u_phase = phase[i][0,mode_index]
                self.u_phase.append((phase[i][0,mode_index] - u_phase) % 360)
                self.v_phase.append((phase[i][1,mode_index] - u_phase) % 360)
                self.w_phase.append((phase[i][2,mode_index] - u_phase) % 360)
                
                self.p_phase.append((phase[i][3,mode_index] - u_phase) % 360)
                self.q_phase.append((phase[i][4,mode_index] - u_phase) % 360)
                self.r_phase.append((phase[i][5,mode_index] - u_phase) % 360)
                
                self.x_phase.append((phase[i][6,mode_index] - u_phase) % 360)
                self.y_phase.append((phase[i][7,mode_index] - u_phase) % 360)
                self.z_phase.append((phase[i][8,mode_index] - u_phase) % 360)
                
                self.phi_phase.append((phase[i][9,mode_index] - u_phase) % 360)
                self.theta_phase.append((phase[i][10,mode_index] - u_phase) % 360)
                self.psi_phase.append((phase[i][11,mode_index] - u_phase) % 360)

def parse_matrix(matrix, indices_list):
    new_matrix = np.zeros((len(matrix[:,0]),len(indices_list)))
    
    for i in range(len(indices_list)):
        
        indices = indices_list[i]
        values = []
        
        for row, col in indices:
            try:
                value = matrix[row][col]
                values.append(value)
            except IndexError:
                # Handle case where index is out of bounds
                values.append(None)
                print(f"Index out of bounds: [{row}, {col}]")
        new_matrix[:,i] = np.array(values)
    
    return new_matrix

def parse_eigvec_matrix_new(matrix, indices_list):
    new_matrix = np.zeros((12,8))
    
    for i in range(len(indices_list)):
        
        indices = indices_list[i]
        
        for _, col in indices:
            column = matrix[:,col]
            # try:
            #     column = matrix[:,col]
            #     values.append(value)
            # except IndexError:
            #     # Handle case where index is out of bounds
            #     values.append(None)
            #     print(f"Index out of bounds: [{row}, {col}]")
        new_matrix[:,i] = np.array(column)
    
    return new_matrix

def parse_eigvec_matrix(matrix, indices_list, index):
    new_matrix = np.zeros((12,8))
    
    for i in range(len(indices_list)):
        
        indices = indices_list[i][index]
        

            
        column = matrix[:,indices[1]]
            # try:
            #     column = matrix[:,col]
            #     values.append(value)
            # except IndexError:
            #     # Handle case where index is out of bounds
            #     values.append(None)
            #     print(f"Index out of bounds: [{row}, {col}]")
        new_matrix[:,i] = np.array(column)
    
    return new_matrix
            
def plot_eigvecs(input_array, eigvecs_full, eigvecs_derivs, eigvecs_coords, mode_label = 'none', case_label='_none',  num_tics = None,
                 x_axis_label='Bank Angle ($\phi$), deg', directory='./', plot_amps = True, plot_phase = True, save_figs=False,
                 ylims = [(None,None),(None,None),(None,None),(None,None)], xlim = (None,None)):
    
    
    alpha_coord = 1.0
    alpha_deriv = 1.0

    # color_coord = '0.3'
    # color_deriv = '0.6'
    
    color_coord = 'darkgray'
    color_deriv = 'indianred'
    
    linestyle_1 = '-'
    linestyle_2 = '-'
    linestyle_3 = '-'
    
    marker_size = 15

    figsize = (5.25,5.5)
    
    fig_amp, ax_amp = plt.subplots(2,2,figsize=figsize, sharex=True)
    # plt.locator_params(axis='x', nbins=5)
    if plot_amps == True:
        # plt.figure(figsize=(4,3))
        ax_amp[0,0].scatter(input_array, eigvecs_full.u_amps, marker='x', color='k', s = marker_size, label='$\Delta u$', clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_full.v_amps, marker='o', color='k', s = marker_size, label='$\Delta v$', clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_full.w_amps, marker='s', color='k', s = marker_size, label='$\Delta w$', clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_full.u_amps, color='k', clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_full.v_amps, color='k', clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_full.w_amps, color='k', clip_on=False)
    
        ax_amp[0,0].scatter(input_array, eigvecs_derivs.u_amps, marker='x', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_derivs.v_amps, marker='o', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_derivs.w_amps, marker='s', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_derivs.u_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_derivs.v_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_derivs.w_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
    
        ax_amp[0,0].scatter(input_array, eigvecs_coords.u_amps, marker='x', color=color_coord, s = marker_size, alpha = alpha_deriv, clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_coords.v_amps, marker='o', color=color_coord, s = marker_size, alpha = alpha_deriv, clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_coords.w_amps, marker='s', color=color_coord, s = marker_size, alpha = alpha_deriv, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_coords.u_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_deriv, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_coords.v_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_deriv, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_coords.w_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_deriv, clip_on=False)
    
    
        ax_amp[0,0].grid()
        # ax_amp[0,0].set_xlabel(x_axis_label)
        ax_amp[0,0].set_ylabel('Amplitude')
        ax_amp[0,0].legend(handletextpad=0.00)
        ax_amp[0,0].set_xlim(xlim)
        ax_amp[0,0].set_ylim(ylims[0])
        # ax_amp.tight_layout()
        # ax_amp.show()
        
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_uvw_eigenvectors_amps_' + case_label + '.png', dpi=250)
    
        # ax_amp.figure(figsize=(4,3))
        ax_amp[0,1].scatter(input_array, eigvecs_full.p_amps, marker='x', color='k', s = marker_size, label='$\Delta p$', clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_full.q_amps, marker='o', color='k', s = marker_size, label='$\Delta q$', clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_full.r_amps, marker='s', color='k', s = marker_size, label='$\Delta r$', clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_full.p_amps, color='k', clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_full.q_amps, color='k', clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_full.r_amps, color='k', clip_on=False)
    
        ax_amp[0,1].scatter(input_array, eigvecs_derivs.p_amps, marker='x', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_derivs.q_amps, marker='o', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_derivs.r_amps, marker='s', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_derivs.p_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_derivs.q_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_derivs.r_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
    
        ax_amp[0,1].scatter(input_array, eigvecs_coords.p_amps, marker='x', color=color_coord, s = marker_size, alpha = alpha_coord, clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_coords.q_amps, marker='o', color=color_coord, s = marker_size, alpha = alpha_coord, clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_coords.r_amps, marker='s', color=color_coord, s = marker_size, alpha = alpha_coord, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_coords.p_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_coords.q_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_coords.r_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord, clip_on=False)
    
        ax_amp[0,1].grid()
        # ax_amp[0,1].set_xlabel(x_axis_label)
        # ax_amp[0,1].set_ylabel('Amplitude')
        ax_amp[0,1].legend(handletextpad=0.00)
        ax_amp[0,1].set_xlim(xlim)
        ax_amp[0,1].set_ylim(ylims[1])
        ax_amp[0,1].yaxis.major.formatter.set_powerlimits((-3,3))
        # ax_amp.tight_layout()
        # ax_amp.show()
    
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_pqr_eigenvectors_amps_' + case_label + '.png', dpi=250)
    
        # ax_amp.figure(figsize=(4,3))
        ax_amp[1,0].scatter(input_array, eigvecs_full.x_amps, marker='x', color='k', s = marker_size, label='$\Delta x_c$', clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_full.y_amps, marker='o', color='k', s = marker_size, label='$\Delta y_c$', clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_full.z_amps, marker='s', color='k', s = marker_size, label='$\Delta z_c$', clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_full.x_amps, color='k', clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_full.y_amps, color='k', clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_full.z_amps, color='k', clip_on=False)
    
        ax_amp[1,0].scatter(input_array, eigvecs_derivs.x_amps, marker='x', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_derivs.y_amps, marker='o', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_derivs.z_amps, marker='s', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_derivs.x_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_derivs.y_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_derivs.z_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
    
        ax_amp[1,0].scatter(input_array, eigvecs_coords.x_amps, marker='x', color=color_coord, s = marker_size, alpha = alpha_coord, clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_coords.y_amps, marker='o', color=color_coord, s = marker_size, alpha = alpha_coord, clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_coords.z_amps, marker='s', color=color_coord, s = marker_size, alpha = alpha_coord, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_coords.x_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_coords.y_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_coords.z_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord, clip_on=False)
    
        ax_amp[1,0].grid()
        ax_amp[1,0].set_xlabel(x_axis_label)
        ax_amp[1,0].set_ylabel('Amplitude')
        ax_amp[1,0].legend(handletextpad=0.00)
        ax_amp[1,0].set_xlim(xlim)
        ax_amp[1,0].set_ylim(ylims[2])
        ax_amp[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_xyz_eigenvectors_amps_' + case_label + '.png', dpi=250)
            
        # ax_amp.figure(figsize=(4,3))
        ax_amp[1,1].scatter(input_array, eigvecs_full.phi_amps, marker='x', color='k', s = marker_size, label='$\Delta\phi$', clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_full.theta_amps, marker='o', color='k', s = marker_size, label='$\Delta\\theta$', clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_full.psi_amps, marker='s', color='k', s = marker_size, label='$\Delta\psi$', clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_full.phi_amps, color='k', clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_full.theta_amps, color='k', clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_full.psi_amps, color='k', clip_on=False)
    
        ax_amp[1,1].scatter(input_array, eigvecs_derivs.phi_amps, marker='x', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_derivs.theta_amps, marker='o', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_derivs.psi_amps, marker='s', color=color_deriv, s = marker_size, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_derivs.phi_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_derivs.theta_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_derivs.psi_amps, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv, clip_on=False)
    
        ax_amp[1,1].scatter(input_array, eigvecs_coords.phi_amps, marker='x', color=color_coord, s = marker_size, alpha = alpha_coord, clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_coords.theta_amps, marker='o', color=color_coord, s = marker_size, alpha = alpha_coord, clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_coords.psi_amps, marker='s', color=color_coord, s = marker_size, alpha = alpha_coord, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_coords.phi_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_coords.theta_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_coords.psi_amps, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord, clip_on=False)
    
        ax_amp[1,1].grid()
        ax_amp[1,1].set_xlabel(x_axis_label)
        # ax_amp[1,1].set_ylabel('Amplitude')
        ax_amp[1,1].legend(handletextpad=0.00)
        ax_amp[1,1].set_xlim(xlim)
        ax_amp[1,1].set_ylim(ylims[3])
        ax_amp[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        ax_amp[1,1].yaxis.major.formatter.set_powerlimits((-3,3))
        
        plt.tight_layout()
        plt.show()
    
        if save_figs == True:
            plt.savefig(directory + mode_label + '_eigenvectors_amps_' + case_label + '.pdf')

    
    if plot_phase == True:
            
        fig_phase, ax_phase = plt.subplots(2,2,figsize=figsize)
        
        # plt.figure(figsize=(4,3))
        ax_phase[0,0].scatter(input_array, eigvecs_full.u_phase, marker='x', color='k', s = marker_size, label='$\Delta u$')
        ax_phase[0,0].scatter(input_array, eigvecs_full.v_phase, marker='o', color='k', s = marker_size, label='$\Delta v$')
        ax_phase[0,0].scatter(input_array, eigvecs_full.w_phase, marker='s', color='k', s = marker_size, label='$\Delta w$')
        ax_phase[0,0].plot(input_array, eigvecs_full.u_phase, color='k')
        ax_phase[0,0].plot(input_array, eigvecs_full.v_phase, color='k')
        ax_phase[0,0].plot(input_array, eigvecs_full.w_phase, color='k')
    
        ax_phase[0,0].scatter(input_array, eigvecs_derivs.u_phase, marker='x', color=color_deriv, s = marker_size, alpha= alpha_deriv)
        ax_phase[0,0].scatter(input_array, eigvecs_derivs.v_phase, marker='o', color=color_deriv, s = marker_size, alpha= alpha_deriv)
        ax_phase[0,0].scatter(input_array, eigvecs_derivs.w_phase, marker='s', color=color_deriv, s = marker_size, alpha= alpha_deriv)
        ax_phase[0,0].plot(input_array, eigvecs_derivs.u_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv)
        ax_phase[0,0].plot(input_array, eigvecs_derivs.v_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv)
        ax_phase[0,0].plot(input_array, eigvecs_derivs.w_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv)
    
        ax_phase[0,0].scatter(input_array, eigvecs_coords.u_phase, marker='x', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[0,0].scatter(input_array, eigvecs_coords.v_phase, marker='o', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[0,0].scatter(input_array, eigvecs_coords.w_phase, marker='s', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[0,0].plot(input_array, eigvecs_coords.u_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
        ax_phase[0,0].plot(input_array, eigvecs_coords.v_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
        ax_phase[0,0].plot(input_array, eigvecs_coords.w_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
    
    
        ax_phase[0,0].grid()
        ax_phase[0,0].set_xlabel(x_axis_label)
        ax_phase[0,0].set_ylabel('Phase [deg]')
        ax_phase[0,0].legend(handletextpad=0.00)
        # plt.tight_layout()
        # plt.show()
        
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_uvw_eigenvectors_phase_' + case_label + '.png', dpi=250)
    
        # plt.figure(figsize=(4,3))
        ax_phase[0,1].scatter(input_array, eigvecs_full.p_phase, marker='x', color='k', s = marker_size, label='$\Delta p$')
        ax_phase[0,1].scatter(input_array, eigvecs_full.q_phase, marker='o', color='k', s = marker_size, label='$\Delta q$')
        ax_phase[0,1].scatter(input_array, eigvecs_full.r_phase, marker='s', color='k', s = marker_size, label='$\Delta r$')
        ax_phase[0,1].plot(input_array, eigvecs_full.p_phase, color='k')
        ax_phase[0,1].plot(input_array, eigvecs_full.q_phase, color='k')
        ax_phase[0,1].plot(input_array, eigvecs_full.r_phase, color='k')
    
        ax_phase[0,1].scatter(input_array, eigvecs_derivs.p_phase, marker='x', color=color_deriv, s = marker_size, alpha= alpha_deriv)
        ax_phase[0,1].scatter(input_array, eigvecs_derivs.q_phase, marker='o', color=color_deriv, s = marker_size, alpha= alpha_deriv)
        ax_phase[0,1].scatter(input_array, eigvecs_derivs.r_phase, marker='s', color=color_deriv, s = marker_size, alpha= alpha_deriv)
        ax_phase[0,1].plot(input_array, eigvecs_derivs.p_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv)
        ax_phase[0,1].plot(input_array, eigvecs_derivs.q_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv)
        ax_phase[0,1].plot(input_array, eigvecs_derivs.r_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv)
    
        ax_phase[0,1].scatter(input_array, eigvecs_coords.p_phase, marker='x', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[0,1].scatter(input_array, eigvecs_coords.q_phase, marker='o', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[0,1].scatter(input_array, eigvecs_coords.r_phase, marker='s', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[0,1].plot(input_array, eigvecs_coords.p_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
        ax_phase[0,1].plot(input_array, eigvecs_coords.q_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
        ax_phase[0,1].plot(input_array, eigvecs_coords.r_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
    
        ax_phase[0,1].grid()
        ax_phase[0,1].set_xlabel(x_axis_label)
        ax_phase[0,1].set_ylabel('Phase [deg]')
        ax_phase[0,1].legend(handletextpad=0.00)
        # plt.tight_layout()
        # plt.show()
    
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_pqr_eigenvectors_phase_' + case_label + '.png', dpi=250)
    
        # plt.figure(figsize=(4,3))
        ax_phase[1,0].scatter(input_array, eigvecs_full.x_phase, marker='x', color='k', s = marker_size, label='$\Delta x_c$')
        ax_phase[1,0].scatter(input_array, eigvecs_full.y_phase, marker='o', color='k', s = marker_size, label='$\Delta y_c$')
        ax_phase[1,0].scatter(input_array, eigvecs_full.z_phase, marker='s', color='k', s = marker_size, label='$\Delta z_c$')
        ax_phase[1,0].plot(input_array, eigvecs_full.x_phase, color='k')
        ax_phase[1,0].plot(input_array, eigvecs_full.y_phase, color='k')
        ax_phase[1,0].plot(input_array, eigvecs_full.z_phase, color='k')
    
        ax_phase[1,0].scatter(input_array, eigvecs_derivs.x_phase, marker='x', color=color_deriv, s = marker_size, alpha= alpha_deriv)
        ax_phase[1,0].scatter(input_array, eigvecs_derivs.y_phase, marker='o', color=color_deriv, s = marker_size, alpha= alpha_deriv)
        ax_phase[1,0].scatter(input_array, eigvecs_derivs.z_phase, marker='s', color=color_deriv, s = marker_size, alpha= alpha_deriv)
        ax_phase[1,0].plot(input_array, eigvecs_derivs.x_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv)
        ax_phase[1,0].plot(input_array, eigvecs_derivs.y_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv)
        ax_phase[1,0].plot(input_array, eigvecs_derivs.z_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_deriv)
    
        ax_phase[1,0].scatter(input_array, eigvecs_coords.x_phase, marker='x', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[1,0].scatter(input_array, eigvecs_coords.y_phase, marker='o', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[1,0].scatter(input_array, eigvecs_coords.z_phase, marker='s', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[1,0].plot(input_array, eigvecs_coords.x_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
        ax_phase[1,0].plot(input_array, eigvecs_coords.y_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
        ax_phase[1,0].plot(input_array, eigvecs_coords.z_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
    
        ax_phase[1,0].grid()
        ax_phase[1,0].set_xlabel(x_axis_label)
        ax_phase[1,0].set_ylabel('Phase [deg]')
        ax_phase[1,0].legend(handletextpad=0.00)
        # plt.tight_layout()
        # plt.show()
    
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_xyz_eigenvectors_phase_' + case_label + '.png', dpi=250)
            
        # plt.figure(figsize=(4,3))
        ax_phase[1,1].scatter(input_array, eigvecs_full.phi_phase, marker='x', color='k', s = marker_size, label='$\Delta\phi$')
        ax_phase[1,1].scatter(input_array, eigvecs_full.theta_phase, marker='o', color='k', s = marker_size, label='$\Delta\\theta$')
        ax_phase[1,1].scatter(input_array, eigvecs_full.psi_phase, marker='s', color='k', s = marker_size, label='$\Delta\psi$')
        ax_phase[1,1].plot(input_array, eigvecs_full.phi_phase, color='k')
        ax_phase[1,1].plot(input_array, eigvecs_full.theta_phase, color='k')
        ax_phase[1,1].plot(input_array, eigvecs_full.psi_phase, color='k')
    
        ax_phase[1,1].scatter(input_array, eigvecs_derivs.phi_phase, marker='x', color=color_deriv, s = marker_size, alpha= alpha_coord)
        ax_phase[1,1].scatter(input_array, eigvecs_derivs.theta_phase, marker='o', color=color_deriv, s = marker_size, alpha= alpha_coord)
        ax_phase[1,1].scatter(input_array, eigvecs_derivs.psi_phase, marker='s', color=color_deriv, s = marker_size, alpha= alpha_coord)
        ax_phase[1,1].plot(input_array, eigvecs_derivs.phi_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_coord)
        ax_phase[1,1].plot(input_array, eigvecs_derivs.theta_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_coord)
        ax_phase[1,1].plot(input_array, eigvecs_derivs.psi_phase, color=color_deriv, linestyle=linestyle_2, alpha= alpha_coord)
    
        ax_phase[1,1].scatter(input_array, eigvecs_coords.phi_phase, marker='x', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[1,1].scatter(input_array, eigvecs_coords.theta_phase, marker='o', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[1,1].scatter(input_array, eigvecs_coords.psi_phase, marker='s', color=color_coord, s = marker_size, alpha = alpha_coord)
        ax_phase[1,1].plot(input_array, eigvecs_coords.phi_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
        ax_phase[1,1].plot(input_array, eigvecs_coords.theta_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
        ax_phase[1,1].plot(input_array, eigvecs_coords.psi_phase, color=color_coord, linestyle=linestyle_3, alpha = alpha_coord)
    
        ax_phase[1,1].grid()
        ax_phase[1,1].set_xlabel(x_axis_label)
        ax_phase[1,1].set_ylabel('Phase [deg]')
        ax_phase[1,1].legend(handletextpad=0.00)
        plt.tight_layout()
        plt.show()
    
        if save_figs == True:
            plt.savefig(directory + mode_label + '_eigenvectors_phase_' + case_label + '.pdf')
        
def plot_single_eigvec_amplitudes(input_array, eigvecs_full, color = 'k', mode_label = 'none',
                                  case_label='_none', x_axis_label='Bank Angle, $\phi$ [deg]',
                                 save_figs = False, directory = './', xlim = None,
                                 ylims = [(None,None),(None,None),(None,None),(None,None)]):
    
    alpha_coord = 1.0
    alpha_deriv = 1.0

    color_coord = 'r'
    color_deriv = 'b'

    figsize = (5.25,5.5)
    fig_amps, ax_amps = plt.subplots(2,2,figsize=figsize,sharex=True)
    # plt.locator_params(axis='x', nbins=5)

    # plt.figure(figsize=(4,3))
    ax_amps[0,0].scatter(input_array, eigvecs_full.u_amps, marker='x', color=color, s = 30, label='$\Delta u$', clip_on=False)
    ax_amps[0,0].scatter(input_array, eigvecs_full.v_amps, marker='o', color=color, s = 30, label='$\Delta v$', clip_on=False)
    ax_amps[0,0].scatter(input_array, eigvecs_full.w_amps, marker='s', color=color, s = 30, label='$\Delta w$', clip_on=False)
    ax_amps[0,0].plot(input_array, eigvecs_full.u_amps, color=color, clip_on=False)
    ax_amps[0,0].plot(input_array, eigvecs_full.v_amps, color=color, clip_on=False)
    ax_amps[0,0].plot(input_array, eigvecs_full.w_amps, color=color, clip_on=False)

    ax_amps[0,0].grid(visible=True)
    # ax_amps[0,0].set_xlabel(x_axis_label)
    ax_amps[0,0].set_ylabel('Amplitude')
    ax_amps[0,0].legend(handletextpad=0.00)
    ax_amps[0,0].set_xlim(xlim)
    ax_amps[0,0].set_ylim(ylims[0])
    # plt.tight_layout()
    # plt.show()
    
    # if save_figs == True:
    #     plt.savefig(directory + mode_label + '_uvw_eigenvectors_' + case_label + '.png', dpi=250)

    # plt.figure(figsize=(4,3))
    ax_amps[0,1].scatter(input_array, eigvecs_full.p_amps, marker='x', color=color, s = 30, label='$\Delta p$', clip_on=False)
    ax_amps[0,1].scatter(input_array, eigvecs_full.q_amps, marker='o', color=color, s = 30, label='$\Delta q$', clip_on=False)
    ax_amps[0,1].scatter(input_array, eigvecs_full.r_amps, marker='s', color=color, s = 30, label='$\Delta r$', clip_on=False)
    ax_amps[0,1].plot(input_array, eigvecs_full.p_amps, color=color, clip_on=False)
    ax_amps[0,1].plot(input_array, eigvecs_full.q_amps, color=color, clip_on=False)
    ax_amps[0,1].plot(input_array, eigvecs_full.r_amps, color=color, clip_on=False)

    ax_amps[0,1].grid(visible=True)
    # ax_amps[0,1].set_xlabel(x_axis_label)
    # ax_amps[0,1].set_ylabel('Amplitude')
    ax_amps[0,1].legend(handletextpad=0.00)
    ax_amps[0,1].set_xlim(xlim)
    ax_amps[0,1].set_ylim(ylims[1])
    # plt.tight_layout()
    # plt.show()

    # if save_figs == True:
    #     plt.savefig(directory + mode_label + '_pqr_eigenvectors_' + case_label + '.png', dpi=250)

    # plt.figure(figsize=(4,3))
    ax_amps[1,0].scatter(input_array, eigvecs_full.x_amps, marker='x', color=color, s = 30, label='$\Delta x_c$', clip_on=False)
    ax_amps[1,0].scatter(input_array, eigvecs_full.y_amps, marker='o', color=color, s = 30, label='$\Delta y_c$', clip_on=False)
    ax_amps[1,0].scatter(input_array, eigvecs_full.z_amps, marker='s', color=color, s = 30, label='$\Delta z_c$', clip_on=False)
    ax_amps[1,0].plot(input_array, eigvecs_full.x_amps, color=color, clip_on=False)
    ax_amps[1,0].plot(input_array, eigvecs_full.y_amps, color=color, clip_on=False)
    ax_amps[1,0].plot(input_array, eigvecs_full.z_amps, color=color, clip_on=False)


    ax_amps[1,0].grid(visible=True)
    ax_amps[1,0].set_xlabel(x_axis_label)
    ax_amps[1,0].set_ylabel('Amplitude')
    ax_amps[1,0].legend(handletextpad=0.00)
    ax_amps[1,0].set_xlim(xlim)
    ax_amps[1,0].set_ylim(ylims[2])
    # ax_amps[1,0].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax_amps[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
    # plt.tight_layout()
    # plt.show()

    # if save_figs == True:
    #     plt.savefig(directory + mode_label + '_xyz_eigenvectors_' + case_label + '.png', dpi=250)
        
    # plt.figure(figsize=(4,3))
    ax_amps[1,1].scatter(input_array, eigvecs_full.phi_amps, marker='x', color=color, s = 30, label='$\Delta\phi$', clip_on=False)
    ax_amps[1,1].scatter(input_array, eigvecs_full.theta_amps, marker='o', color=color, s = 30, label='$\Delta\\theta$', clip_on=False)
    ax_amps[1,1].scatter(input_array, eigvecs_full.psi_amps, marker='s', color=color, s = 30, label='$\Delta\psi$', clip_on=False)
    ax_amps[1,1].plot(input_array, eigvecs_full.phi_amps, color=color, clip_on=False)
    ax_amps[1,1].plot(input_array, eigvecs_full.theta_amps, color=color, clip_on=False)
    ax_amps[1,1].plot(input_array, eigvecs_full.psi_amps, color=color, clip_on=False)


    ax_amps[1,1].grid(visible=True)
    ax_amps[1,1].set_xlabel(x_axis_label)
    # ax_amps[1,1].set_ylabel('Amplitude')
    ax_amps[1,1].legend(handletextpad=0.00)
    ax_amps[1,1].set_xlim(xlim)
    ax_amps[1,1].set_ylim(ylims[3])
    # ax_amps[1,1].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax_amps[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
    plt.tight_layout()
    plt.show()

    if save_figs == True:
        plt.savefig(directory + mode_label + '_eigenvectors_amps_' + case_label + '.pdf')

def plot_polar_eigvecs(amps,phase,state_vars = ['u','v','w']):
   
    fig, ax = plt.subplots(1,1,figsize=(3,3), subplot_kw={"projection" : "polar"})
    
    lbl_params = dict(ha="center",va="center",size=10.0)
    bbox_dict = dict(facecolor="w",linewidth=0,alpha=1.0, boxstyle="Square, pad=0.05")
    
    label_dict = {'u':0,'v':1,'w':2,'p':3,'q':4,'r':5,'x':6,'y':7,'z':8,'phi':9,'theta':10,'psi':11}

    labels = ['$u$','$v$','$w$','$p$','$q$','$r$','$x$','$y$','$z$','$\phi$','$\\theta$','$\psi$']

    indices = [label_dict[key] for key in state_vars if key in label_dict]
    markers = ['o','s','8','^','<','>','x','1','4','P','p','*']
    
    index_max = np.argmax(amps[indices])
    
    for i in range(len(indices)):
        index = indices[i]
        if index == index_max:
            # ax.plot([0., phase[index]],[0., amps[index]],lw=2.0,color='r',marker=markers[index],markevery=[-1],label=labels[index])
            ax.plot([0., phase[index]],[0., amps[index]],lw=2.0,color='r',label=labels[index])
            # ax.text(phase[index],amps[index]*1.0,'T', bbox = bbox_dict, **lbl_params)
        else:
            ax.plot([0., phase[index]],[0., amps[indices][index_max]],lw=0.9,color='k', linestyle = '--')
            # ax.plot([0., phase[index]],[0., amps[index]],lw=2.0,color='r',marker=markers[index],markevery=[-1],label=labels[index])
            ax.plot([0., phase[index]],[0., amps[index]],lw=2.0,color='r',label=labels[index])
            # ax.text(phase[index],amps[index]*1.0,'T', bbox = bbox_dict, **lbl_params)
            # ax.text(phase[index],amps[index]*1.0,'T', bbox = bbox_dict, **lbl_params)
        
    # ax.grid(True)
    ax.set_rlabel_position(60.0)
    # ax.set_xlabel("Phase Angle [deg]")
    # plt.legend()
    angle = np.deg2rad(67.5)
    # ax.legend(loc="lower left", bbox_to_anchor=(.5 + np.cos(angle)/2, .5 + np.sin(angle)/2))
    plt.tight_layout()
    # ax[0,0].text(np.deg2rad(90.0),0.3*max(A_p),"Amplitude",rotation=60.0)
    plt.show()

def plot_eigenvalues_imag_LWD_RWD(eig_real, eig_imag, xlims, ylims, sim_results = None, num_tics = None,
                          lin_sim_results = None, legend_on = False,
                          save_fig = False, filename = 'temp', directory = './'):
    
    '''FUNCTION FOR COMPARING LWD AND RWD BOOMERANG SOLUTIONS'''
    
    fig_size = (3.25, 3.25)
    label1 = 'RWD'
    label2 = 'LWD'
    
    marker = matplotlib.markers.MarkerStyle(marker='s', fillstyle='none')
    
    axes_linewidth = 0.6
    
    xaxis_label = 'Real Component'
    yaxis_label = 'Imaginary Component'
    size = 5

    plt.figure(figsize=fig_size)
    plt.get_current_fig_manager().set_window_title(filename)

    for i in range(len(eig_real[0][:])):
        
        if i == len(eig_real[0][:]) - 1:
            labela = label1
            labelb = label2
        else:
            labela = ''
            labelb = '' 

        "RWD"
        plt.scatter(eig_real[0][i], eig_imag[0][i], marker=marker, color='k', s = size,label=labela, clip_on=False)
        
        "LWD"
        if len(eig_real) == 2:
            plt.scatter(eig_real[1][i], eig_imag[1][i], marker=marker, color='r', s = size,label=labelb, clip_on=False)

        size += 3.5

    plt.axhline(y=0, color='k',linewidth=axes_linewidth)
    plt.grid(visible=True)

    plt.xlabel(xaxis_label)
    plt.ylabel(yaxis_label)
    plt.ylim(ylims)
    plt.xlim(xlims)
    ax = plt.gca()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
    
    if legend_on:
        plt.legend()
        
    plt.tight_layout()
    plt.show()
    
    if save_fig == True:
        plt.savefig(directory + filename + '.pdf')
        plt.savefig(directory + filename + '.png',dpi=400)

def plot_eigenvalues_real_LWD_RWD(x_values, eig_real, xlims, ylims, num_tics = None, legend_on = False, x_axis_label = 'x-axis label',
                                 save_fig = False, filename = 'temp', directory = './'):
    fig_size = (3.25, 3.25)

    label1 = 'RWD'
    label2 = 'LWD'
    
    marker = matplotlib.markers.MarkerStyle(marker='s', fillstyle='none')

    axes_linewidth = 0.6
    y_axis_label = 'Real Component'
    size = 30
    plt.figure(figsize=fig_size)
    '''SPIRAL'''

    plt.scatter(x_values, eig_real[0], marker=marker, color='k', s = size,label=label1, clip_on=False)
    if len(eig_real) == 2:
        plt.scatter(x_values, eig_real[1], marker=marker, color='r', s = size,label=label2, clip_on=False)

    plt.axhline(y=0, color='k',linewidth=axes_linewidth)
    plt.grid(visible=True)
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    if legend_on:
        plt.legend()
    ax = plt.gca()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
    plt.ylim(ylims)
    plt.xlim(xlims)
    plt.tight_layout()
    plt.show()
    
    if save_fig == True:
        plt.savefig(directory + filename + '.pdf')
        plt.savefig(directory + filename + '.png',dpi=400)

def plot_eigenvalues_imag(eig_real, eig_imag, xlims, ylims, sim_results = None,
                          lin_sim_results = None, legend_on = False, num_tics = None,
                          save_fig = False, filename = 'temp', directory = './'):
    fig_size = (3.25, 3.25)
    # label1 = 'Eq. (26)'
    # label2 = 'Classical C.S.'
    # label3 = 'Symmetric Aircraft'
    # label4 = 'Nonlinear SIM'
    # label5 = 'Linear SIM'
    
    label1 = 'Eq. (27)'
    label2 = 'COORD'
    label3 = 'SYM'
    label4 = 'NL-SIM'
    label5 = 'L-SIM'
    
    marker = matplotlib.markers.MarkerStyle(marker='s', fillstyle='none')
    DERIVS_marker = matplotlib.markers.MarkerStyle(marker='o', fillstyle='none')
    COORD_marker = matplotlib.markers.MarkerStyle(marker='^', fillstyle='none')
    SIM_marker = matplotlib.markers.MarkerStyle(marker='x', fillstyle='none')
    
    SIM_COLOR = 'indianred'
    LIN_SIM_COLOR = 'royalblue'
    
    
    axes_linewidth = 0.6
    
    xaxis_label = 'Real Component'
    yaxis_label = 'Imaginary Component'
    size = 5

    plt.figure(figsize=fig_size)

    for i in range(len(eig_real[0][:])):
        
        if i == len(eig_real[0][:]) - 1:
            labela = label1
            labelb = label2
            labelc = label3
            labeld = label4
            labele = label5
        else:
            labela = ''
            labelb = '' 
            labelc = ''
            labeld = ''
            labele = ''

        plt.scatter(eig_real[0][i], eig_imag[0][i], marker=marker, color='k', s = size,label=labela, clip_on=False)

        '''Symetric Aircraft Approximations'''
        plt.scatter(eig_real[1][i], eig_imag[1][i], marker=DERIVS_marker, color='k', s = size,label=labelc, clip_on=False)

        '''Coord System Alignment Approximation'''
        plt.scatter(eig_real[2][i], eig_imag[2][i], marker=COORD_marker, color='k', s = size,label=labelb, clip_on=False)

        if sim_results != None:
            plt.scatter(sim_results[0][i],sim_results[1][i], marker=SIM_marker, color=SIM_COLOR, s = size, label=labeld, clip_on=False)
        if lin_sim_results != None:
            plt.scatter(lin_sim_results[0][i],lin_sim_results[1][i], marker=SIM_marker, color=LIN_SIM_COLOR, s = size, label=labele, clip_on=False)

        size += 3.5

    plt.axhline(y=0, color='k',linewidth=axes_linewidth)
    plt.grid(visible=True)

    plt.xlabel(xaxis_label)
    plt.ylabel(yaxis_label)
    plt.ylim(ylims)
    plt.xlim(xlims)
    ax = plt.gca()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
    
    if legend_on:
        plt.legend()
        
    plt.tight_layout()
    plt.show()
    
    if save_fig == True:
        plt.savefig(directory + filename + '.pdf')
    
def plot_eigenvalues_imag_single(eig_real, eig_imag, xlims, ylims, legend_on = False,
                                 save_fig = False, filename = 'temp', directory = './'):
    fig_size = (3.25, 3.25)
    label1 = 'Eq. (26)'
    marker = matplotlib.markers.MarkerStyle(marker='s', fillstyle='none')
    axes_linewidth = 0.6
    
    xaxis_label = 'Real Component'
    yaxis_label = 'Imaginary Component'
    size = 5

    plt.figure(figsize=fig_size)

    for i in range(len(eig_real[:])):
        
        if i == len(eig_real[:]) - 1:
            labela = label1
        else:
            labela = ''

        plt.scatter(eig_real[i], eig_imag[i], marker=marker, color='k', s = size,label=labela, clip_on=False)

        size += 3.5

    plt.axhline(y=0, color='k',linewidth=axes_linewidth)
    plt.grid(visible=True)

    plt.xlabel(xaxis_label)
    plt.ylabel(yaxis_label)
    plt.ylim(ylims)
    plt.xlim(xlims)
    
    ax = plt.gca()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    
    if legend_on:
        plt.legend()
    plt.tight_layout()
    plt.show()
    
    if save_fig == True:
        plt.savefig(directory + filename + '.pdf')
    
def plot_eigenvalues_real(x_values, eig_real, xlims, ylims, legend_on = False, x_axis_label = 'x-axis label', num_tics = None,
                                 save_fig = False, filename = 'temp', directory = './'):
    fig_size = (3.25, 3.25)
    # label1 = 'Eq. (26)'
    # label2 = 'Classical C.S.'
    # label3 = 'Symmetric Aircraft'
    
    label1 = 'Eq. (27)'
    label2 = 'COORD'
    label3 = 'SYM'
    
    marker = matplotlib.markers.MarkerStyle(marker='s', fillstyle='none')
    DERIVS_marker = matplotlib.markers.MarkerStyle(marker='o', fillstyle='none')
    COORD_marker = matplotlib.markers.MarkerStyle(marker='^', fillstyle='none')

    axes_linewidth = 0.6
    y_axis_label = 'Real Component'
    size = 30
    plt.figure(figsize=fig_size)
    '''SPIRAL'''

    plt.scatter(x_values, eig_real[0], marker=marker, color='k', s = size,label=label1, clip_on=False)
    plt.scatter(x_values, eig_real[1], marker=DERIVS_marker, color='k', s = size,label=label3, clip_on=False)
    plt.scatter(x_values, eig_real[2], marker=COORD_marker, color='k', s = size,label=label2, clip_on=False)

    plt.axhline(y=0, color='k',linewidth=axes_linewidth)
    plt.grid(visible=True)
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    plt.xlim(xlims)
    plt.ylim(ylims)

    ax = plt.gca()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
    plt.tight_layout()
    plt.show()
    
    if legend_on:
        plt.legend()
    
    if save_fig == True:
        plt.savefig(directory + filename + '.pdf')
    
def plot_eigenvalues_real_single(x_values, eig_real, xlims, ylims, legend_on = False, x_axis_label = 'x-axis label',
                                 save_fig = False, filename = 'temp', directory = './'):
    fig_size = (3.25, 3.25)
    label1 = 'Eq. (26)'

    marker = matplotlib.markers.MarkerStyle(marker='s', fillstyle='none')
    axes_linewidth = 0.6
    y_axis_label = 'Real Component'
    size = 30

    plt.figure(figsize=fig_size)
    '''SPIRAL'''

    plt.scatter(x_values, eig_real, marker=marker, color='k', s = size,label=label1, clip_on=False)

    plt.axhline(y=0, color='k',linewidth=axes_linewidth)
    plt.grid(visible=True)
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    if legend_on:
        plt.legend()
    ax = plt.gca()
    # ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    plt.tight_layout()
    plt.show()
    
    if save_fig == True:
        plt.savefig(directory + filename + '.pdf')

'''BIRE VS F16 COMPARISONS'''

def plot_eigenvalues_imag_double(eig_real, eig_imag, xlims, ylims, sim_results = None,
                          lin_sim_results = None, legend_on = False,
                          save_fig = False, filename = 'temp', directory = './'):
    fig_size = (3.25, 3.25)
    label1 = 'Baseline'
    label2 = 'BIRE'
    
    marker = matplotlib.markers.MarkerStyle(marker='s', fillstyle='none')
    BIRE_marker = matplotlib.markers.MarkerStyle(marker='o', fillstyle='none')
    
    axes_linewidth = 0.6
    
    xaxis_label = 'Real Component'
    yaxis_label = 'Imaginary Component'
    size = 5

    plt.figure(figsize=fig_size)

    for i in range(len(eig_real[0][:])):
        
        if i == len(eig_real[0][:]) - 1:
            labela = label1
            labelb = label2
        else:
            labela = ''
            labelb = '' 

        plt.scatter(eig_real[0][i], eig_imag[0][i], marker=marker, color='k', s = size,label=labela, clip_on=False)
        
        '''BIRE'''
        plt.scatter(eig_real[1][i], eig_imag[1][i], marker=BIRE_marker, color='k', s = size,label=labelb, clip_on=False)

        size += 3.5

    plt.axhline(y=0, color='k',linewidth=axes_linewidth)
    plt.grid(visible=True)

    plt.xlabel(xaxis_label)
    plt.ylabel(yaxis_label)
    plt.ylim(ylims)
    plt.xlim(xlims)
    ax = plt.gca()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    
    if legend_on:
        plt.legend()
        
    plt.tight_layout()
    plt.show()
    
    if save_fig == True:
        plt.savefig(directory + filename + '.pdf')

def plot_eigenvalues_real_double(x_values, eig_real, xlims, ylims, legend_on = False, x_axis_label = 'x-axis label',
                                 save_fig = False, filename = 'temp', directory = './'):
    fig_size = (3.25, 3.25)
    label1 = 'Baseline'
    label2 = 'BIRE'
    
    marker = matplotlib.markers.MarkerStyle(marker='s', fillstyle='none')
    BIRE_marker = matplotlib.markers.MarkerStyle(marker='o', fillstyle='none')
    
    axes_linewidth = 0.6
    y_axis_label = 'Real Component'
    size = 30
    plt.figure(figsize=fig_size)
    '''SPIRAL'''

    plt.scatter(x_values, eig_real[0], marker=marker, color='k', s = size,label=label1, clip_on=False)
    plt.scatter(x_values, eig_real[1], marker=BIRE_marker, color='k', s = size,label=label2, clip_on=False)

    plt.axhline(y=0, color='k',linewidth=axes_linewidth)
    plt.grid(visible=True)
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    if legend_on:
        plt.legend()
    plt.ylim(ylims)
    plt.xlim(xlims)
    ax = plt.gca()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    # ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
    # ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    # ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.tight_layout()
    plt.show()
    
    if save_fig == True:
        plt.savefig(directory + filename + '.pdf')

'''BASELINE VS BIRE EIGENVECTOR PLOTS'''
    
def plot_eigvecs_double(input_array, eigvecs_base, eigvecs_BIRE, xlim = (0,60), mode_label = 'none', case_label='_none',
                       x_axis_label='Bank Angle ($\phi$), deg', directory='./', plot_amps = True, plot_phase = True, save_figs=False,
                       ylims = [(None,None),(None,None),(None,None),(None,None)]):
    
    alpha_base = 1.0
    alpha_BIRE = 1.0
    
    linestyle_1 = '-'
    linestyle_2 = '-'
    linestyle_3 = '-'

    color_BIRE = 'k'
    color_base = 'darkgray'

    figsize = (5.25,5.5)
    
    fig_amp, ax_amp = plt.subplots(2,2,figsize=figsize, sharex=True)
    # plt.locator_params(axis='x', nbins=5)
    if plot_amps == True:
        # plt.figure(figsize=(4,3))
        ax_amp[0,0].scatter(input_array, eigvecs_base.u_amps, marker='x', color=color_base, s = 30, clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_base.v_amps, marker='o', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_base.w_amps, marker='s', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_base.u_amps, color=color_base, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_base.v_amps, color=color_base, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_base.w_amps, color=color_base, clip_on=False)
    
        ax_amp[0,0].scatter(input_array, eigvecs_BIRE.u_amps, marker='x', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta u$', clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_BIRE.v_amps, marker='o', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta v$', clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_BIRE.w_amps, marker='s', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta w$', clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_BIRE.u_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_BIRE.v_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_BIRE.w_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
    

        ax_amp[0,0].grid()
        # ax_amp[0,0].set_xlabel(x_axis_label)
        ax_amp[0,0].set_ylabel('Amplitude')
        ax_amp[0,0].legend(handletextpad=0.00)
        # ax_amp[0,0].legend(handletextpad=0.00, mode = "expand", ncol = 3)
        ax_amp[0,0].set_xlim(xlim)
        ax_amp[0,0].set_ylim(ylims[0])
        # ax_amp[0,0].sharex(ax_amp[1,0])
        # ax_amp[0,0].set_xticks([])
        # ax_amp.tight_layout()
        # ax_amp.show()
        
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_uvw_eigenvectors_amps_' + case_label + '.png', dpi=250)
    
        # ax_amp.figure(figsize=(4,3))
        ax_amp[0,1].scatter(input_array, eigvecs_base.p_amps, marker='x', color=color_base, s = 30, clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_base.q_amps, marker='o', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_base.r_amps, marker='s', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_base.p_amps, color=color_base, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_base.q_amps, color=color_base, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_base.r_amps, color=color_base, clip_on=False)
    
        ax_amp[0,1].scatter(input_array, eigvecs_BIRE.p_amps, marker='x', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta p$', clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_BIRE.q_amps, marker='o', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta q$', clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_BIRE.r_amps, marker='s', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta r$', clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_BIRE.p_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_BIRE.q_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_BIRE.r_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
    
    
        ax_amp[0,1].grid()
        # ax_amp[0,1].set_xlabel(x_axis_label)
        # ax_amp[0,1].set_ylabel('Amplitude')
        ax_amp[0,1].legend(handletextpad=0.00)
        ax_amp[0,1].set_xlim(xlim)
        ax_amp[0,1].set_ylim(ylims[1])
        # ax_amp[0,1].sharex(ax_amp[1,1])
        # ax_amp[0,1].set_xticks([])
        # plt.locator_params(axis='x', nbins=5)
        # ax_amp.tight_layout()
        # ax_amp.show()
    
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_pqr_eigenvectors_amps_' + case_label + '.png', dpi=250)
    
        # ax_amp.figure(figsize=(4,3))
        ax_amp[1,0].scatter(input_array, eigvecs_base.x_amps, marker='x', color=color_base, s = 30, clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_base.y_amps, marker='o', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_base.z_amps, marker='s', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_base.x_amps, color=color_base, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_base.y_amps, color=color_base, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_base.z_amps, color=color_base, clip_on=False)
    
        ax_amp[1,0].scatter(input_array, eigvecs_BIRE.x_amps, marker='x', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta x_c$', clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_BIRE.y_amps, marker='o', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta y_c$', clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_BIRE.z_amps, marker='s', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta z_c$', clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_BIRE.x_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_BIRE.y_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_BIRE.z_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
    
        ax_amp[1,0].grid()
        ax_amp[1,0].set_xlabel(x_axis_label)
        ax_amp[1,0].set_ylabel('Amplitude')
        ax_amp[1,0].legend(handletextpad=0.00)
        ax_amp[1,0].set_xlim(xlim)
        ax_amp[1,0].set_ylim(ylims[2])
        # ax_amp[1,0].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax_amp[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
        # ax_amp.tight_layout()
        # ax_amp.show()
    
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_xyz_eigenvectors_amps_' + case_label + '.png', dpi=250)
            
        # ax_amp.figure(figsize=(4,3))
        ax_amp[1,1].scatter(input_array, eigvecs_base.phi_amps, marker='x', color=color_base, s = 30, clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_base.theta_amps, marker='o', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_base.psi_amps, marker='s', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_base.phi_amps, color=color_base, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_base.theta_amps, color=color_base, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_base.psi_amps, color=color_base, clip_on=False)
    
        ax_amp[1,1].scatter(input_array, eigvecs_BIRE.phi_amps, marker='x', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta\phi$', clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_BIRE.theta_amps, marker='o', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta\\theta$', clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_BIRE.psi_amps, marker='s', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta\psi$', clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_BIRE.phi_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_BIRE.theta_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_BIRE.psi_amps, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
    
        ax_amp[1,1].grid()
        ax_amp[1,1].set_xlabel(x_axis_label)
        # ax_amp[1,1].set_ylabel('Amplitude')
        ax_amp[1,1].legend(handletextpad=0.00)
        ax_amp[1,1].set_xlim(xlim)
        ax_amp[1,1].set_ylim(ylims[3])
        # ax_amp[1,1].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax_amp[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
        # plt.locator_params(axis='x', nbins=5)
        plt.tight_layout()
        plt.show()
    
        if save_figs == True:
            plt.savefig(directory + mode_label + '_eigenvectors_amps_' + case_label + '.pdf')
    
    
    if plot_phase == True:
        fig_phase, ax_phase = plt.subplots(2,2,figsize=figsize)        
        # plt.figure(figsize=(4,3))
        ax_phase[0,0].scatter(input_array, eigvecs_base.u_phase, marker='x', color=color_base, s = 30, clip_on=False)
        ax_phase[0,0].scatter(input_array, eigvecs_base.v_phase, marker='o', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_phase[0,0].scatter(input_array, eigvecs_base.w_phase, marker='s', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_phase[0,0].plot(input_array, eigvecs_base.u_phase, color=color_base, clip_on=False)
        ax_phase[0,0].plot(input_array, eigvecs_base.v_phase, color=color_base, clip_on=False)
        ax_phase[0,0].plot(input_array, eigvecs_base.w_phase, color=color_base, clip_on=False)
    
        ax_phase[0,0].scatter(input_array, eigvecs_BIRE.u_phase, marker='x', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta u$', clip_on=False)
        ax_phase[0,0].scatter(input_array, eigvecs_BIRE.v_phase, marker='o', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta v$', clip_on=False)
        ax_phase[0,0].scatter(input_array, eigvecs_BIRE.w_phase, marker='s', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta w$', clip_on=False)
        ax_phase[0,0].plot(input_array, eigvecs_BIRE.u_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_phase[0,0].plot(input_array, eigvecs_BIRE.v_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_phase[0,0].plot(input_array, eigvecs_BIRE.w_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
    
        ax_phase[0,0].grid()
        ax_phase[0,0].set_xlabel(x_axis_label)
        ax_phase[0,0].set_ylabel('Phase [deg]')
        ax_phase[0,0].legend(handletextpad=0.00)
        ax_phase[0,0].set_xlim(xlim)
        # ax_phase[0,0].set_ylim(ylims[0])
        # plt.tight_layout()
        # plt.show()
        
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_uvw_eigenvectors_phase_' + case_label + '.png', dpi=250)
    
        # plt.figure(figsize=(4,3))
        ax_phase[0,1].scatter(input_array, eigvecs_base.p_phase, marker='x', color=color_base, s = 30, clip_on=False)
        ax_phase[0,1].scatter(input_array, eigvecs_base.q_phase, marker='o', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_phase[0,1].scatter(input_array, eigvecs_base.r_phase, marker='s', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_phase[0,1].plot(input_array, eigvecs_base.p_phase, color=color_base, clip_on=False)
        ax_phase[0,1].plot(input_array, eigvecs_base.q_phase, color=color_base, clip_on=False)
        ax_phase[0,1].plot(input_array, eigvecs_base.r_phase, color=color_base, clip_on=False)
    
        ax_phase[0,1].scatter(input_array, eigvecs_BIRE.p_phase, marker='x', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta p$', clip_on=False)
        ax_phase[0,1].scatter(input_array, eigvecs_BIRE.q_phase, marker='o', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta q$', clip_on=False)
        ax_phase[0,1].scatter(input_array, eigvecs_BIRE.r_phase, marker='s', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta r$', clip_on=False)
        ax_phase[0,1].plot(input_array, eigvecs_BIRE.p_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_phase[0,1].plot(input_array, eigvecs_BIRE.q_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_phase[0,1].plot(input_array, eigvecs_BIRE.r_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)

        ax_phase[0,1].grid()
        ax_phase[0,1].set_xlabel(x_axis_label)
        ax_phase[0,1].set_ylabel('Phase [deg]')
        ax_phase[0,1].legend(handletextpad=0.00)
        ax_phase[0,1].set_xlim(xlim)
        # ax_phase[0,1].set_ylim(ylims[1])
        # plt.tight_layout()
        # plt.show()
    
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_pqr_eigenvectors_phase_' + case_label + '.png', dpi=250)
    
        # plt.figure(figsize=(4,3))
        ax_phase[1,0].scatter(input_array, eigvecs_base.x_phase, marker='x', color=color_base, s = 30, clip_on=False)
        ax_phase[1,0].scatter(input_array, eigvecs_base.y_phase, marker='o', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_phase[1,0].scatter(input_array, eigvecs_base.z_phase, marker='s', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_phase[1,0].plot(input_array, eigvecs_base.x_phase, color=color_base, clip_on=False)
        ax_phase[1,0].plot(input_array, eigvecs_base.y_phase, color=color_base, clip_on=False)
        ax_phase[1,0].plot(input_array, eigvecs_base.z_phase, color=color_base, clip_on=False)
    
        ax_phase[1,0].scatter(input_array, eigvecs_BIRE.x_phase, marker='x', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta x_c$', clip_on=False)
        ax_phase[1,0].scatter(input_array, eigvecs_BIRE.y_phase, marker='o', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta y_c$', clip_on=False)
        ax_phase[1,0].scatter(input_array, eigvecs_BIRE.z_phase, marker='s', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta z_c$', clip_on=False)
        ax_phase[1,0].plot(input_array, eigvecs_BIRE.x_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_phase[1,0].plot(input_array, eigvecs_BIRE.y_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_phase[1,0].plot(input_array, eigvecs_BIRE.z_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
    
        ax_phase[1,0].grid()
        ax_phase[1,0].set_xlabel(x_axis_label)
        ax_phase[1,0].set_ylabel('Phase [deg]')
        ax_phase[1,0].legend(handletextpad=0.00)
        ax_phase[1,0].set_xlim(xlim)
        # ax_phase[1,0].set_ylim(ylims[2])
        # plt.tight_layout()
        # plt.show()
    
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_xyz_eigenvectors_phase_' + case_label + '.png', dpi=250)
            
        # plt.figure(figsize=(4,3))
        ax_phase[1,1].scatter(input_array, eigvecs_base.phi_phase, marker='x', color=color_base, s = 30, clip_on=False)
        ax_phase[1,1].scatter(input_array, eigvecs_base.theta_phase, marker='o', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_phase[1,1].scatter(input_array, eigvecs_base.psi_phase, marker='s', color=color_base, facecolors='none', s = 30, clip_on=False)
        ax_phase[1,1].plot(input_array, eigvecs_base.phi_phase, color=color_base, clip_on=False)
        ax_phase[1,1].plot(input_array, eigvecs_base.theta_phase, color=color_base, clip_on=False)
        ax_phase[1,1].plot(input_array, eigvecs_base.psi_phase, color=color_base, clip_on=False)
    
        ax_phase[1,1].scatter(input_array, eigvecs_BIRE.phi_phase, marker='x', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta\phi$', clip_on=False)
        ax_phase[1,1].scatter(input_array, eigvecs_BIRE.theta_phase, marker='o', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta\\theta$', clip_on=False)
        ax_phase[1,1].scatter(input_array, eigvecs_BIRE.psi_phase, marker='s', color=color_BIRE, s = 30, alpha= alpha_BIRE, label='$\Delta\psi$', clip_on=False)
        ax_phase[1,1].plot(input_array, eigvecs_BIRE.phi_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_phase[1,1].plot(input_array, eigvecs_BIRE.theta_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
        ax_phase[1,1].plot(input_array, eigvecs_BIRE.psi_phase, color=color_BIRE, linestyle=linestyle_2, alpha= alpha_BIRE, clip_on=False)
    
        ax_phase[1,1].grid()
        ax_phase[1,1].set_xlabel(x_axis_label)
        ax_phase[1,1].set_ylabel('Phase [deg]')
        ax_phase[1,1].legend(handletextpad=0.00)
        ax_phase[1,1].set_xlim(xlim)
        # ax_phase[1,1].set_ylim(ylims[3])
        plt.tight_layout()
        plt.show()
    
        if save_figs == True:
            plt.savefig(directory + mode_label + '_eigenvectors_phase_' + case_label + '.pdf')

'''PLOT LWD RWD EIGENVECTOR SOLUTIONS FOR THE BOOMERANG'''
        
def plot_eigvecs_LWD_RWD(input_array, eigvecs_RWD, eigvecs_LWD = None, mode_label = 'none', case_label='_none', num_tics = None,
                       x_axis_label='Bank Angle, $\phi$, deg', directory='./', plot_amps = True, plot_phase = True, save_figs=False,
                       xlim = (0,60), ylims = [(None,None),(None,None),(None,None),(None,None)]):
        
    alpha_RWD = 1.0
    alpha_LWD = 1.0

    color_RWD = 'k'
    color_LWD = 'r'
    
    linestyle_1 = '-'
    linestyle_2 = '-'
    linestyle_3 = '-'

    figsize = (5.25,5.5)
    marker_size = 30
    
    # figsize = (3.25,3.5)
    # marker_size = 15
    fig_amp, ax_amp = plt.subplots(2,2,figsize=figsize,sharex=True)
    
    if plot_amps == True:
        # plt.figure(figsize=(4,3))
        ax_amp[0,0].scatter(input_array, eigvecs_RWD.u_amps, marker='x', color=color_RWD, s = marker_size, label='$\Delta u$', clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_RWD.v_amps, marker='o', color=color_RWD, s = marker_size, label='$\Delta v$', clip_on=False)
        ax_amp[0,0].scatter(input_array, eigvecs_RWD.w_amps, marker='s', color=color_RWD, s = marker_size, label='$\Delta w$', clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_RWD.u_amps, color=color_RWD, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_RWD.v_amps, color=color_RWD, clip_on=False)
        ax_amp[0,0].plot(input_array, eigvecs_RWD.w_amps, color=color_RWD, clip_on=False)  
            
        if eigvecs_LWD != None:
            ax_amp[0,0].scatter(input_array, eigvecs_LWD.u_amps, marker='x', color=color_LWD, s = marker_size, clip_on=False)
            ax_amp[0,0].scatter(input_array, eigvecs_LWD.v_amps, marker='o', color=color_LWD, facecolors='none', s = marker_size, clip_on=False)
            ax_amp[0,0].scatter(input_array, eigvecs_LWD.w_amps, marker='s', color=color_LWD, facecolors='none', s = marker_size, clip_on=False)
            ax_amp[0,0].plot(input_array, eigvecs_LWD.u_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)
            ax_amp[0,0].plot(input_array, eigvecs_LWD.v_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)
            ax_amp[0,0].plot(input_array, eigvecs_LWD.w_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)    
    
        ax_amp[0,0].grid()
        # ax_amp[0,0].set_xlabel(x_axis_label)
        ax_amp[0,0].set_ylabel('Amplitude')
        ax_amp[0,0].legend(handletextpad=0.00)
        ax_amp[0,0].set_xlim(xlim)
        ax_amp[0,0].set_ylim(ylims[0])
        # ax_amp.tight_layout()
        # ax_amp.show()
        
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_uvw_eigenvectors_amps_' + case_label + '.png', dpi=250)
    
        # ax_amp.figure(figsize=(4,3))
        ax_amp[0,1].scatter(input_array, eigvecs_RWD.p_amps, marker='x', color=color_RWD, s = marker_size, label='$\Delta p$', clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_RWD.q_amps, marker='o', color=color_RWD, s = marker_size, label='$\Delta q$', clip_on=False)
        ax_amp[0,1].scatter(input_array, eigvecs_RWD.r_amps, marker='s', color=color_RWD, s = marker_size, label='$\Delta r$', clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_RWD.p_amps, color=color_RWD, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_RWD.q_amps, color=color_RWD, clip_on=False)
        ax_amp[0,1].plot(input_array, eigvecs_RWD.r_amps, color=color_RWD, clip_on=False)
            
        if eigvecs_LWD != None:
            ax_amp[0,1].scatter(input_array, eigvecs_LWD.p_amps, marker='x', color=color_LWD, s = marker_size, clip_on=False)
            ax_amp[0,1].scatter(input_array, eigvecs_LWD.q_amps, marker='o', color=color_LWD, facecolors='none', s = marker_size, clip_on=False)
            ax_amp[0,1].scatter(input_array, eigvecs_LWD.r_amps, marker='s', color=color_LWD, facecolors='none', s = marker_size, clip_on=False)
            ax_amp[0,1].plot(input_array, eigvecs_LWD.p_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)
            ax_amp[0,1].plot(input_array, eigvecs_LWD.q_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)
            ax_amp[0,1].plot(input_array, eigvecs_LWD.r_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)    
    
    
        ax_amp[0,1].grid()
        # ax_amp[0,1].set_xlabel(x_axis_label)
        # ax_amp[0,1].set_ylabel('Amplitude')
        ax_amp[0,1].legend(handletextpad=0.00)
        ax_amp[0,1].set_xlim(xlim)
        ax_amp[0,1].set_ylim(ylims[1])
        ax_amp[0,1].yaxis.major.formatter.set_powerlimits((-3,3))
        # ax_amp.tight_layout()
        # ax_amp.show()
    
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_pqr_eigenvectors_amps_' + case_label + '.png', dpi=250)
    
        # ax_amp.figure(figsize=(4,3))
        ax_amp[1,0].scatter(input_array, eigvecs_RWD.x_amps, marker='x', color=color_RWD, s = marker_size, label='$\Delta x_c$', clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_RWD.y_amps, marker='o', color=color_RWD, s = marker_size, label='$\Delta y_c$', clip_on=False)
        ax_amp[1,0].scatter(input_array, eigvecs_RWD.z_amps, marker='s', color=color_RWD, s = marker_size, label='$\Delta z_c$', clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_RWD.x_amps, color=color_RWD, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_RWD.y_amps, color=color_RWD, clip_on=False)
        ax_amp[1,0].plot(input_array, eigvecs_RWD.z_amps, color=color_RWD, clip_on=False)
            
        if eigvecs_LWD != None:
            ax_amp[1,0].scatter(input_array, eigvecs_LWD.x_amps, marker='x', color=color_LWD, s = marker_size, clip_on=False)
            ax_amp[1,0].scatter(input_array, eigvecs_LWD.y_amps, marker='o', color=color_LWD, facecolors='none', s = marker_size, clip_on=False)
            ax_amp[1,0].scatter(input_array, eigvecs_LWD.z_amps, marker='s', color=color_LWD, facecolors='none', s = marker_size, clip_on=False)
            ax_amp[1,0].plot(input_array, eigvecs_LWD.x_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)
            ax_amp[1,0].plot(input_array, eigvecs_LWD.y_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)
            ax_amp[1,0].plot(input_array, eigvecs_LWD.z_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)    
    
        ax_amp[1,0].grid()
        ax_amp[1,0].set_xlabel(x_axis_label)
        ax_amp[1,0].set_ylabel('Amplitude')
        ax_amp[1,0].legend(handletextpad=0.00)
        ax_amp[1,0].set_xlim(xlim)
        ax_amp[1,0].set_ylim(ylims[2])
        ax_amp[1,0].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        # ax_amp.tight_layout()
        # ax_amp.show()
    
        # if save_figs == True:
        #     plt.savefig(directory + mode_label + '_xyz_eigenvectors_amps_' + case_label + '.png', dpi=250)
            
        # ax_amp.figure(figsize=(4,3))
        ax_amp[1,1].scatter(input_array, eigvecs_RWD.phi_amps, marker='x', color=color_RWD, s = marker_size, label='$\Delta\phi$', clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_RWD.theta_amps, marker='o', color=color_RWD, s = marker_size, label='$\Delta\\theta$', clip_on=False)
        ax_amp[1,1].scatter(input_array, eigvecs_RWD.psi_amps, marker='s', color=color_RWD,  s = marker_size, label='$\Delta\psi$', clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_RWD.phi_amps, color=color_RWD, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_RWD.theta_amps, color=color_RWD, clip_on=False)
        ax_amp[1,1].plot(input_array, eigvecs_RWD.psi_amps, color=color_RWD, clip_on=False)
            
        if eigvecs_LWD != None:
            ax_amp[1,1].scatter(input_array, eigvecs_LWD.phi_amps, marker='x', color=color_LWD, s = marker_size, clip_on=False)
            ax_amp[1,1].scatter(input_array, eigvecs_LWD.theta_amps, marker='o', color=color_LWD, facecolors='none', s = marker_size, clip_on=False)
            ax_amp[1,1].scatter(input_array, eigvecs_LWD.psi_amps, marker='s', color=color_LWD, facecolors='none', s = marker_size, clip_on=False)
            ax_amp[1,1].plot(input_array, eigvecs_LWD.phi_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)
            ax_amp[1,1].plot(input_array, eigvecs_LWD.theta_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)
            ax_amp[1,1].plot(input_array, eigvecs_LWD.psi_amps, color=color_LWD, linestyle=linestyle_3, clip_on=False)    
    
        ax_amp[1,1].grid()
        ax_amp[1,1].set_xlabel(x_axis_label)
        # ax_amp[1,1].set_ylabel('Amplitude')
        ax_amp[1,1].legend(handletextpad=0.00)
        ax_amp[1,1].set_xlim(xlim)
        ax_amp[1,1].set_ylim(ylims[3])
        ax_amp[1,1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=num_tics))
        ax_amp[1,1].yaxis.major.formatter.set_powerlimits((-3,3))
        
        plt.tight_layout()
        plt.show()
    
        if save_figs == True:
            plt.savefig(directory + mode_label + '_eigenvectors_amps_' + case_label + '.pdf')
            plt.savefig(directory + mode_label + '_eigenvectors_amps_' + case_label + '.png',dpi=400)
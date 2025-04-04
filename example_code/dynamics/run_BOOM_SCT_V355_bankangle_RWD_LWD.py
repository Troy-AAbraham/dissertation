from dynamics_analysis import dynamicAnalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# from matplotlib.ticker import FormatStrFormatter

from handle_eigenvectors import *

# options for saving result figures
case_label = 'BOOM_SCT_V355_RWD_LWD'
directory = './results/'
save_figs = True

# initialize lists for results
real_list_RWD = []
imag_list_RWD = []
ev_amps_RWD = []
ev_phase_RWD = []

alpha_RWD_array = []
beta_RWD_array = []
tau_RWD_array = []
da_RWD_array = []
de_RWD_array = []
dr_RWD_array = []
p0_RWD_array = []
q0_RWD_array = []
r0_RWD_array = []

real_list_LWD = []
imag_list_LWD = []
ev_amps_LWD = []
ev_phase_LWD = []

alpha_LWD_array = []
beta_LWD_array = []
tau_LWD_array = []
da_LWD_array = []
de_LWD_array = []
dr_LWD_array = []
p0_LWD_array = []
q0_LWD_array = []
r0_LWD_array = []

derivatives_RWD = derivative_storage(nondim_derivs=False)
derivatives_LWD = derivative_storage(nondim_derivs=False)

'''------- Set Flight Condition -------'''

V = 355 #ft/s
gamma = np.deg2rad(0.0) #rad
H = 24000. #ft
cg_shift = [0., 0., 0.] #ft

SHSS = False
COMP = False
STALL = False

# bank angle range
bank_angle_array = np.linspace(0.0,60.0,11)

'''RUN RWD SCT RANGE'''

for i in range(len(bank_angle_array)):
    phi = np.deg2rad(bank_angle_array[i]) #rad

    case = dynamicAnalysis(path='./', write_output = False, output_filename = 'eig_vals_boomerang.txt',
                            shss=SHSS, compressible=COMP, stall=STALL,
                            coords_approx= False, derivs_approx=False, cg_shift=cg_shift)
    case.solve_equilibrium_state(V, H, gamma, phi) 
    case.solve_derivatives()    
    case.solve_dynamics_system()
    
    # append results
    
    derivatives_RWD.append_derivatives(case)
    
    tau, alpha, beta, da, de, dr = case.eq_inputs
    alpha_RWD_array.append(alpha)
    beta_RWD_array.append(beta)
    da_RWD_array.append(da)
    de_RWD_array.append(de)
    dr_RWD_array.append(dr)
    tau_RWD_array.append(tau)
    
    p0_RWD_array.append(case.p0)
    q0_RWD_array.append(case.q0)
    r0_RWD_array.append(case.r0)
        
    real_list_RWD.append(case.eigreal)
    imag_list_RWD.append(case.eigimag)
    
    ev_amps_RWD.append(case.amps)
    ev_phase_RWD.append(case.phase)    

# convert to numpy arrays
real_list_RWD = np.asarray(real_list_RWD)
imag_list_RWD = np.asarray(imag_list_RWD)

'''RUN LWD SCT RANGE'''

for i in range(len(bank_angle_array)):
    
    phi = -np.deg2rad(bank_angle_array[i]) #rad

    case = dynamicAnalysis(path='./', write_output = False, output_filename = 'eig_vals_boomerang.txt',
                            shss=SHSS, compressible=COMP, stall=STALL,
                            coords_approx= False, derivs_approx=False, cg_shift=cg_shift)
    case.solve_equilibrium_state(V, H, gamma, phi) 
    case.solve_derivatives()    
    case.solve_dynamics_system()
    
    derivatives_LWD.append_derivatives(case)
    
    tau, alpha, beta, da, de, dr = case.eq_inputs
    alpha_LWD_array.append(alpha)
    beta_LWD_array.append(beta)
    da_LWD_array.append(da)
    de_LWD_array.append(de)
    dr_LWD_array.append(dr)
    tau_LWD_array.append(tau)
    
    p0_LWD_array.append(case.p0)
    q0_LWD_array.append(case.q0)
    r0_LWD_array.append(case.r0)
        
    real_list_LWD.append(case.eigreal)
    imag_list_LWD.append(case.eigimag)
    
    ev_amps_LWD.append(case.amps)
    ev_phase_LWD.append(case.phase)
    
real_list_LWD = np.asarray(real_list_LWD)
imag_list_LWD = np.asarray(imag_list_LWD)


'''------ Plot Results -------'''

# set indices for modes for parsing results
LP_ind = [6,7]
SP_ind = [10,11]
DR_ind = [8,9]
SR_ind = 4
RL_ind = 5

# some organization for plotting

LP_eigvecs_RWD = parse_eigenvectors(ev_amps_RWD, ev_phase_RWD, LP_ind[0])
SP_eigvecs_RWD = parse_eigenvectors(ev_amps_RWD, ev_phase_RWD, SP_ind[0])
DR_eigvecs_RWD = parse_eigenvectors(ev_amps_RWD, ev_phase_RWD, DR_ind[0])
RL_eigvecs_RWD = parse_eigenvectors(ev_amps_RWD, ev_phase_RWD, RL_ind)
SR_eigvecs_RWD = parse_eigenvectors(ev_amps_RWD, ev_phase_RWD, SR_ind)

LP_eigvecs_LWD = parse_eigenvectors(ev_amps_LWD, ev_phase_LWD, LP_ind[0])
SP_eigvecs_LWD = parse_eigenvectors(ev_amps_LWD, ev_phase_LWD, SP_ind[0])
DR_eigvecs_LWD = parse_eigenvectors(ev_amps_LWD, ev_phase_LWD, DR_ind[0])
RL_eigvecs_LWD = parse_eigenvectors(ev_amps_LWD, ev_phase_LWD, RL_ind)
SR_eigvecs_LWD = parse_eigenvectors(ev_amps_LWD, ev_phase_LWD, SR_ind)


LP_real_eig_vals = [real_list_RWD[:,LP_ind[0]], real_list_LWD[:,LP_ind[0]]]
LP_imag_eig_vals = [imag_list_RWD[:,LP_ind[0]], imag_list_LWD[:,LP_ind[0]]]
SP_real_eig_vals = [real_list_RWD[:,SP_ind[0]], real_list_LWD[:,SP_ind[0]]]
SP_imag_eig_vals = [imag_list_RWD[:,SP_ind[0]], imag_list_LWD[:,SP_ind[0]]]
DR_real_eig_vals = [real_list_RWD[:,DR_ind[0]], real_list_LWD[:,DR_ind[0]]]
DR_imag_eig_vals = [imag_list_RWD[:,DR_ind[0]], imag_list_LWD[:,DR_ind[0]]]
RL_eig_vals = [real_list_RWD[:,RL_ind], real_list_LWD[:,RL_ind]]
SR_eig_vals = [real_list_RWD[:,SR_ind], real_list_LWD[:,SR_ind]]

# set limits in eigenvalue figures

xlim_LP = (-0.15,0.075)
ylim_LP = (0.1,0.26)
xlim_SP = (-2.5,-1.5)
ylim_SP = (7.4,8.8)
xlim_DR = (-0.45,-0.41)
ylim_DR = (3.48,3.6)
ylim_RL = (-3.6,-2.2)
ylim_SR = (-0.2, 0.1)

xlim = (0,60)

f_size = 8
plt.rcParams.update({'font.size': f_size})

plot_eigenvalues_imag_LWD_RWD(eig_real=LP_real_eig_vals, eig_imag = LP_imag_eig_vals,
                            xlims = xlim_LP, ylims = ylim_LP, legend_on=True, num_tics = 4,
                            filename = 'LP_eigenvalue_compare_' + case_label, directory = directory + 'eigenvalues/', save_fig = save_figs)

plot_eigenvalues_imag_LWD_RWD(eig_real=SP_real_eig_vals, eig_imag = SP_imag_eig_vals,
                            xlims = xlim_SP, ylims = ylim_SP, legend_on=False, num_tics = 5,
                            filename = 'SP_eigenvalue_compare_' + case_label, directory = directory + 'eigenvalues/', save_fig = save_figs)

plot_eigenvalues_imag_LWD_RWD(eig_real=DR_real_eig_vals, eig_imag = DR_imag_eig_vals,
                            xlims = xlim_DR, ylims = ylim_DR, legend_on=False, num_tics = 5,
                            filename = 'DR_eigenvalue_compare_' + case_label, directory = directory + 'eigenvalues/', save_fig = save_figs)

plot_eigenvalues_real_LWD_RWD(x_values = bank_angle_array, eig_real = RL_eig_vals, num_tics = 5,
                            xlims = xlim, ylims = ylim_RL, legend_on=True, x_axis_label = 'Bank Angle ($\phi$), deg',
                            filename = 'RL_eigenvalue_compare_' + case_label, directory = directory + 'eigenvalues/', save_fig = save_figs)

plot_eigenvalues_real_LWD_RWD(x_values = bank_angle_array, eig_real = SR_eig_vals, num_tics = 5,
                            xlims = xlim, ylims = ylim_SR, x_axis_label = 'Bank Angle ($\phi$), deg',
                            filename = 'SR_eigenvalue_compare_' + case_label, directory = directory + 'eigenvalues/', save_fig = save_figs)

# set limits in eigenvector figures

LP_ylims = [(0.0,0.1), (0.0,2e-4), (0.0,1.0), (0.0,1e-3)]
SP_ylims = [(0.0,1.0), (0.0,0.025), (0.0,0.02), (0.0,3e-3)]
DR_ylims = [(0.0,1.0), (0.0,0.01), (0.0,0.015), (0.0,3e-3)]
RL_ylims = [(0.0,1.0), (0.0,0.3), (0.0,1.0), (0.0,0.1)]
SR_ylims = [(0.0,0.1), (0.0,2.5e-4), (0.0,1.0), (0.0,1e-3)]

plot_eigvecs_LWD_RWD(input_array = bank_angle_array, eigvecs_RWD = LP_eigvecs_RWD,
                        eigvecs_LWD = LP_eigvecs_LWD,
                        mode_label = 'LP',
                        case_label = case_label + '_compare',
                        directory=directory + 'eigenvectors/',
                        save_figs = save_figs,
                        ylims = LP_ylims,
                        xlim = xlim,
                        x_axis_label='Bank Angle ($\phi$), deg',
                        plot_phase=False,
                        num_tics = 5)

plot_eigvecs_LWD_RWD(input_array = bank_angle_array, eigvecs_RWD = SP_eigvecs_RWD,
                        eigvecs_LWD = SP_eigvecs_LWD,
                        mode_label = 'SP',
                        case_label = case_label + '_compare',
                        directory=directory + 'eigenvectors/',
                        save_figs = save_figs,
                        ylims = SP_ylims,
                        xlim = xlim,
                        x_axis_label='Bank Angle ($\phi$), deg',
                        plot_phase=False,
                        num_tics = 5)

plot_eigvecs_LWD_RWD(input_array = bank_angle_array, eigvecs_RWD = DR_eigvecs_RWD,
                        eigvecs_LWD = DR_eigvecs_LWD,
                        mode_label = 'DR',
                        case_label = case_label + '_compare',
                        directory=directory + 'eigenvectors/',
                        save_figs = save_figs,
                        ylims = DR_ylims,
                        xlim = xlim,
                        x_axis_label='Bank Angle ($\phi$), deg',
                        plot_phase=False,
                        num_tics = 5)

plot_eigvecs_LWD_RWD(input_array = bank_angle_array, eigvecs_RWD = RL_eigvecs_RWD,
                        eigvecs_LWD = RL_eigvecs_LWD,
                        mode_label = 'RL',
                        case_label = case_label + '_compare',
                        directory=directory + 'eigenvectors/',
                        save_figs = save_figs,
                        ylims = RL_ylims,
                        xlim = xlim,
                        x_axis_label='Bank Angle ($\phi$), deg',
                        plot_phase=False,
                        num_tics = 5)

plot_eigvecs_LWD_RWD(input_array = bank_angle_array, eigvecs_RWD = SR_eigvecs_RWD,
                        eigvecs_LWD = SR_eigvecs_LWD,
                        mode_label = 'SR',
                        case_label = case_label + '_compare',
                        directory=directory + 'eigenvectors/',
                        save_figs = save_figs,
                        ylims = SR_ylims,
                        xlim = xlim,
                        x_axis_label='Bank Angle ($\phi$), deg',
                        plot_phase=False,
                        num_tics = 5)

# plot body-fixed derivatives

derivatives_RWD.plot_body_fixed_derivatives(x_array = bank_angle_array, xlims = (0,60), baseline_derivatives = None, derivatives_2 = derivatives_LWD,
                                        x_label = 'Bank Angle ($\phi$), deg', directory = directory, case_label = case_label, save_figs = save_figs,
                                        ylims = [[(-10,30), (-20,50), (-120,40), (-10000,2000)],
                                                                      [(-16,0), (-200,200), (-600,200), (-10000,2000)],
                                                                      [(-150,50), (-500,100), (-100,200), (-4000,1000)]])

# trim solution figures

f_size = 8
plt.rcParams.update({'font.size': f_size})
figsize = (2.5,2.5)

axes_linewidth = 0.6
plt.figure(figsize=figsize)
plt.plot(bank_angle_array, 180.0*np.array(da_RWD_array)/np.pi, color='k', linestyle = '-', label = '$\delta_a$', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(de_RWD_array)/np.pi, color='k', linestyle = '-.', label = '$\delta_e$', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(dr_RWD_array)/np.pi, color='k', linestyle = ':', label = '$\delta_r$', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(da_LWD_array)/np.pi, color='r', linestyle = '-', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(de_LWD_array)/np.pi, color='r', linestyle = '-.', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(dr_LWD_array)/np.pi, color='r', linestyle = ':', clip_on=False)
plt.xlim(xlim)
plt.ylim(-25,0)
# plt.axhline(y=0, color='k',linewidth=axes_linewidth)
plt.grid(visible=True)
plt.xlabel('Bank Angle ($\phi$), deg')
plt.ylabel('Control Surface Deflection, deg')
ax = plt.gca()
ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
plt.legend()
plt.tight_layout()
plt.show()

if save_figs == True:
    plt.savefig(directory + 'control_deflections_' + case_label + '.pdf') 


plt.figure(figsize=figsize)
plt.plot(bank_angle_array, tau_RWD_array, color='k', linestyle = '-', label = '$\\tau$', clip_on=False)
plt.plot(bank_angle_array, tau_LWD_array, color='r', linestyle = '-', clip_on=False)
plt.xlim(xlim)
plt.ylim(0.5,1.0)
# plt.axhline(y=0, color='k',linewidth=axes_linewidth)
plt.grid(visible=True)
plt.xlabel('Bank Angle ($\phi$), deg')
plt.ylabel('Throttle Setting ($\\tau$)')
ax = plt.gca()
ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
# plt.legend()
plt.tight_layout()
plt.show()

if save_figs == True:
    plt.savefig(directory + 'throttle_setting_' + case_label + '.pdf') 
    
axes_linewidth = 0.6
plt.figure(figsize=figsize)
plt.plot(bank_angle_array, 180.0*np.array(alpha_RWD_array)/np.pi, color='k', linestyle = '-', label = '$\\alpha$', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(beta_RWD_array)/np.pi, color='k', linestyle = '-.', label = '$\\beta$', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(alpha_LWD_array)/np.pi, color='r', linestyle = '-', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(beta_LWD_array)/np.pi, color='r', linestyle = '-.', clip_on=False)

plt.grid(visible=True)
plt.xlabel('Bank Angle ($\phi$), deg')
plt.ylabel('Aerodynamic Angle, deg')
plt.xlim(xlim)
plt.ylim(-5,15)
ax = plt.gca()
ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
plt.legend()
plt.tight_layout()
plt.show()

if save_figs == True:
    plt.savefig(directory + 'aerodynamic_angles_' + case_label + '.pdf') 

plt.figure(figsize=figsize)
plt.plot(bank_angle_array, 180.0*np.array(p0_RWD_array)/np.pi, color='k', linestyle = '-', label = '$p$', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(q0_RWD_array)/np.pi, color='k', linestyle = '-.', label = '$q$', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(r0_RWD_array)/np.pi, color='k', linestyle = ':', label = '$r$', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(p0_LWD_array)/np.pi, color='r', linestyle = '-', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(q0_LWD_array)/np.pi, color='r', linestyle = '-.', clip_on=False)
plt.plot(bank_angle_array, 180.0*np.array(r0_LWD_array)/np.pi, color='r', linestyle = ':', clip_on=False)
plt.xlim(xlim)
plt.ylim(-6,8)

plt.grid(visible=True)
ax = plt.gca()
ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(numticks=5))
plt.xlabel('Bank Angle ($\phi$), deg')
plt.ylabel('Rotation Rate, deg/s')
plt.legend()
plt.tight_layout()
plt.show()

if save_figs == True:
    plt.savefig(directory + 'rotation_rates_' + case_label + '.pdf')
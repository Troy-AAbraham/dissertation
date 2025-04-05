import numpy as np

def body_2_wind_vector(F_b,alpha,beta):
    
    '''Converts a body fixed force vector into a vector in the wind coordinate frame
    
    Parameters
    -----------
    
    F_b: array_like
        body-fixed vector
        
    alpha: float
        angle of attack
        
    beta: float
        sideslip angle
        
    Returns
    -----------
    
    F_w: array-like
        wind-coordinate force vector
    '''
    
    bdy_2_wind_mat = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)],
                               [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(alpha)*np.sin(beta)],
                               [-np.sin(alpha), 0.0, np.cos(alpha)]])
    
    # also changes sign to produce positive lift and drag values
    F_w = np.matmul(bdy_2_wind_mat,F_b)*np.array([-1.0, 1.0, -1.0])
    
    return F_w

def comp_2_body_mat(phi_c,theta_c,psi_c):
    '''
    Creates the direction cosine matrix for the conversion between the 
    component coordinate system and the aircraft coordinate system.
    
    Uses the form that is component to aircraft conversion.
    
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

def comp_velocity(V_body,omega_body,p_body,dir_cos):
    '''
    Calculate the local velocity at a component as a result of aircraft
    velocity, rotational velocity, and component position.
    
    Parameters
    -----------
    
    V_body: array_like
        vehicle body-fixed velocity vector
        
    omega_body: array_like
        vehicle body-fixed roational velocity vector
        
    p_body: array_like
        body-fixed position of the component relative to vehicle CG
        
    Returns
    -----------
    V_comp: array-like
        Velocity vector in the components coordinate system
        
    V_mag: float
        Velocity magnitude at the component
        
    u_comp: array-like
        Unit velocity vector in the components coordinate system
        
    '''
    # Eq. 3.130.3 in Hunsaker, Component contributions to the complete vehicle
    # V = [R]^-1 (V - P_c X w)
    V_comp = np.matmul(np.linalg.inv(dir_cos),(V_body - np.cross(p_body,omega_body)))
    V_mag = np.sqrt(V_comp[0]*V_comp[0] + V_comp[1]*V_comp[1] + V_comp[2]*V_comp[2])
    u_comp = V_comp/V_mag
    
    return V_comp, V_mag, u_comp

def vec_comp_2_body(vector_comp,dir_cosine):
    
    '''
    Converts vector in a components coordinate system to the aircrafts body fixed
    coordinate system.
    
    Parameters
    -----------
    vector: array_like
        component coordinate vector
        
    dir_cosine: array-like (3x3 matrix)
        direction cosine matrix from component to body coordinates
        
    Returns
    -----------
    
    vector_body: array-like
        aircraft body-fixed vector
    '''
    
    vector_body = np.matmul(dir_cosine,vector_comp)
    
    return vector_body

def solve_Vb_vector(alpha, beta, V):
    ''' body fixed velocity components from angle of attack and sideslip 
    angle 
    
    Parameters
    -----------
    alpha: float
        angle of attack
        
    beta: float
        sideslip angle
        
    V: float
        true airspeed

    Returns
    -----------
    u: float
        body-fixed x translational velocity component
        
    v: float
        body-fixed y translational velocity component
        
    w: float
        body-fixed z translational velocity component
    
    '''
    u = V*np.cos(alpha)*np.cos(beta)
    v = V*np.sin(beta)
    w = V*np.sin(alpha)*np.cos(beta)
    return u, v, w

if __name__ == "__main__":

    # Example 3.5.2
    
    V = 200
    
    phi = np.deg2rad(0.0)    
    theta = np.deg2rad(10.0)
    psi = np.deg2rad(5.0)
    
    P_c = np.array([0.1, 5.0, 1.5])
    
    H_c = 2.5
    R_c = 0.5
    
    alpha = 5.0*np.pi/180.
    beta= 5.0*np.pi/180.
    omega = np.array([30.0, 0.0, 0.0])*np.pi/180.
    
    rho = 0.00237689
    mu = 3.737025865723e-7
    
    R = comp_2_body_mat(phi_c = phi, theta_c = theta, psi_c = psi)
    
    Vvec = np.array(solve_Vb_vector(alpha, beta, V))
    
    V_comp, V_mag, u_comp = comp_velocity(V_body = Vvec, omega_body = omega, p_body = P_c, dir_cos = R)
from math import cos, sin, atan, atan2, asin, sqrt, exp, pi
import numpy as np

#-----------------------------------------------------------------------------#
#ATMOSPHERIC PROFILE FUNCTIONS

def gravity_si(H):
   '''Gravitational constant in SI units, m/s^2'''
   return 9.806645*(6356766./(6356766.+H))**2

def gravity_english(H):
   '''Gravitational constant in English units, ft/s^2'''
   return 9.806645*(6356766./(6356766. + H*0.3048))**2/0.3048

def statsi(h):
    '''
    Interpolates atmospheric properties for a standard day. SI units.
    
    Parameters
    -----------
    
    h : float
        altitude in meters

    Returns
    -----------
    h: float
        altitude in meters
        
    z: float
        geopotential altitude in meters
        
    t: float
        absolute temp in Kelvins
        
    p: float
        pressure in N/m^2
        
    d: float
        density in kg/m^3
        
    a: float
        speed of sound in m/s
        
    '''
    Psa = np.zeros(9)
    zsa = [0,11000,20000,32000,47000,52000,61000,79000,90000]
    Tsa = [288.15,216.65,216.65,228.65,270.65,270.65,252.65,180.65,180.65]
    g0 = 9.806645
    R = 287.0528
    Re = 6356766
    gamma = 1.4
    
    Psa[0] = 101325
    z = Re*h/(Re+h)
    for i in range(1,9):
        Lt = -(Tsa[i]-Tsa[i-1])/(zsa[i]-zsa[i-1])
        if Lt == 0:
            if z <= zsa[i]:
                t = Tsa[i-1]
                p = Psa[i-1]*exp(-g0*(z-zsa[i-1])/R/Tsa[i-1])
                d = (p/R)/t
                a = sqrt(gamma*R*t)
                return (h,z,t,p,d,a)
            else:
                Psa[i] = Psa[i-1]*exp(-g0*(zsa[i]-zsa[i-1])/R/Tsa[i-1])
        else:
            ex = (g0/R)/Lt
            if z < zsa[i]:
               t = Tsa[i-1]-Lt*(z-zsa[i-1])
               p = Psa[i-1]*(t/Tsa[i-1])**ex
               d = (p/R)/t
               a = sqrt(gamma*R*t)
               return (h,z,t,p,d,a)
            else:
               Psa[i] = Psa[i-1]*(Tsa[i]/Tsa[i-1])**ex
    t = Tsa[8]
    p = 0.
    d = 0.
    a = sqrt(gamma*R*t)
    
    return (h,z,t,p,d,a)

def statee(h):
    '''
    Converts SI units to english units.
        
    Parameters
    -----------
    
    h : float
        geometric altitude in feet

    Returns
    -----------
    h: float
        geometric altitude in feet
        
    z: float
        geopotential altitude in feet
        
    t: float
        absolute temp in Rakines
        
    p: float
        pressure in lbf/ft^2
        
    d: float
        density in slugs/ft^3
        
    a: float
        speed of sound in ft/s
    '''

    hsi = h*0.3048
    hsi,zsi,tsi,psi,dsi,asi = statsi(hsi)
    z = zsi/0.3048
    t = tsi*1.8
    p = psi*0.020885434304801722
    d = dsi*0.00194032032363104
    a = asi/0.3048
    return (h,z,t,p,d,a)
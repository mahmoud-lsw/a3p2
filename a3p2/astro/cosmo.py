#===========================================================================
# Imports

import math

import numpy as np
import scipy.integrate
import scipy.interpolate

#===========================================================================
# Functions & classes

#---------------------------------------------------------------------------
def cosmo_term1(h=0.7, omega_m=0.3, omega_l=0.7) :
    """
    Returns the cosmological integration terms in units of Gyrs as lambda function.
    """
    return lambda z: 1. / (1.022699E-1 * h) / ((z + 1.) * np.sqrt((z + 1.)**2
                                                                  * (1. + omega_m*z)
                                                                  - z * (2. + z) * omega_l))

#---------------------------------------------------------------------------
def z_to_t_table(zmin=1E-5, zmax=40., nsteps=150, logsteps=True, cosmo=[.7, .3, .7]):
    """
    Calculates cosmic time in Gyrs for a number of redshift points
    and returns them as numpy arrays.
    """
    if logsteps :
        z = np.power(10., np.linspace(math.log10(zmin), math.log10(zmax), nsteps))
    else :
        z = np.linspace(zmin, zmax, nsteps)
    ct1 = cosmo_term1(*cosmo)
    t = np.array([scipy.integrate.quad(ct1, 0., x)[0] for x in z])
    return (z, t);

#===========================================================================

#===========================================================================
# Imports

import math

#===========================================================================
# Functions & classes

#---------------------------------------------------------------------------
def blackbody(T, freq=False, cgssys=False) :
    """
    Returns a black body (Planck's law) as lambda function
    (wavelength/frequency vs spectral radiance).

    T      = Temperature [ Kelvin ]
    freg   = False/(True) - return frequency vs flux (wavelength is standard)
    cgssys = False/(True) - return flux in cgs units (SI is standard)
    """
    from a3p2..constants import si, cgs
    h, c, k = si.h, si.c, si.k
    if cgssys :
        h, c, k = cgs.h, cgs.c, cgs.k
    if freq :
        # SI  : n [ Hz ], F [ J s^-1 m^-2 sr^-1 Hz^-1 ] 
        # CGS : n [ Hz ], F [ erg s^-1 cm^-2 sr^-1 Hz^-1 ]
        return lambda n: 2. * h * n**3. / c**2. / (math.exp(h * n / (k * T)) - 1.)
    else :
        # SI  : l [ m ],  F [ J s^-1 m^-2 sr^-1 m^-1 ] 
        # CGS : l [ cm ], F [ erg s^-1 cm^-2 sr^-1 cm^-1 ]
        return lambda l: 2. * h * c**2. / l**5. / (math.exp(h * c / (l * k * T)) - 1.)

#===========================================================================

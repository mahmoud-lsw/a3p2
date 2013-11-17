#===========================================================================
# Imports

#import sys
import time
import logging
import os
import glob
#import gc
#import gzip
import datetime
#import json

import pyfits
import numpy as np
import scipy.interpolate
import scipy.integrate
import matplotlib.pyplot as plt
#
import a3p2
#import a3p2.astro.cosmo
#import a3p2.tools.config
#import a3p2.ebl.model as eblmodel
#from a3p2.constants import *
import a3p2.ebl.measurements

#===========================================================================

f_lp = lambda x, E0, a, b, phi0: phi0 * (x / E0) ** (-a - b * np.log10(x / E0))

x = 10. ** np.linspace(1., 10., 101)

plt.loglog(x, f_lp(x, 1E3, 0., 5., 1E-6) * x * x)
plt.loglog(x, f_lp(x, 1E7, 0., 5., 1E-14) * x * x)

print scipy.integrate.simps(x=np.log10(x), y=f_lp(x, 1E3, 0., 5., 1E-6) * x * x * np.log(10.))
print scipy.integrate.simps(x=np.log10(x), y=f_lp(x, 1E7, 0., 5., 1E-14) * x * x * np.log(10.))

plt.ylim(1E-30, 1.)

plt.show()
           

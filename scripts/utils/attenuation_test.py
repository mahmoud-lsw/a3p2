#===========================================================================
# Imports

import sys
import time
import logging
import os

import numpy as np
import matplotlib.pyplot as plt

import a3p2
from a3p2.constants import *
from a3p2.ebl.attenuation import *
from a3p2.tools.config import *

#===========================================================================
# Functions

model_file = 'bestfit_BL2006_PopII.model.fits'
ebl = EBLModelMR2012(model_file)
e, z = 1., .2
print e, z
print '1 >>>', calc_ebl_attenuation(z, e, ebl)
print '2 >>>', calc_ebl_attenuation2(z, e, ebl)

nstepsz = [11, 21, 31, 51, 71, 91]
res1, res2, res3, res4 = [], [], [], []
for n in nstepsz :
    res1.append(calc_ebl_attenuation2(z, e, ebl, int_steps_z=21, int_steps_mu=n, int_steps_log10e=501))
    res2.append(calc_ebl_attenuation2(z, e, ebl, int_steps_z=21, int_steps_mu=n, int_steps_log10e=201))
    res3.append(calc_ebl_attenuation2(z, e, ebl, int_steps_z=21, int_steps_mu=n, int_steps_log10e=601))
    res4.append(calc_ebl_attenuation2(z, e, ebl, int_steps_z=21, int_steps_mu=n, int_steps_log10e=1001))    
    #int_steps_z=31, int_steps_mu=21, int_steps_log10e=201

plt.plot(nstepsz, res2, '+', label='201')
plt.plot(nstepsz, res1, 'x', label='501')
plt.plot(nstepsz, res3, 'o', label='601')
plt.plot(nstepsz, res4, 's', label='1001')
plt.legend()
plt.show()

#===========================================================================
#===========================================================================

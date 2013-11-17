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
#import scipy.interpolate
#import scipy.integrate
import matplotlib.pyplot as plt

import a3p2
#import a3p2.astro.cosmo
#import a3p2.tools.config
#import a3p2.ebl.model as eblmodel
from a3p2.constants import *
#import a3p2.ebl.measurements

#===========================================================================
# Functions

#===========================================================================
# Main
filename = 'haardt_madau_2012_UVB.model.fits'

ebl_model_file = a3p2.datapath + '/ebl/models/haardt_madau_2012_UVB.out'
ebl_model_data = np.loadtxt(ebl_model_file)
z = ebl_model_data[0, 1:]
l = ebl_model_data[1:,0] * 1e-4
d = ebl_model_data[1:, 1:] * 3000. / l[:,np.newaxis] * 1E17 * (1. + z[np.newaxis, :]) ** -3
# Sort ascending in frequency (reverse)
d = d[::-1]
l = l[::-1]
nu = si.c / (l * 1E-6)
m = np.diff(nu) > 0.
m = np.hstack((np.array([True]), m))
nu = nu[m]
l = l[m]
d = d[m]

hdulist =  pyfits.HDUList([
    pyfits.PrimaryHDU(),
    pyfits.ImageHDU(data=d,name='EBL'),
    pyfits.new_table(pyfits.ColDefs([
        pyfits.Column(name='LOGNU',format='E', array=np.log10(nu), unit='log10(Hz)')
        ])),
    pyfits.new_table(pyfits.ColDefs([
        pyfits.Column(name='WAVE',format='E', array=l, unit='micrometer')
        ])),
    pyfits.new_table(pyfits.ColDefs([
        pyfits.Column(name='Z',format='E', array=z)
        ]))
    ])
hdulist[-3].header.update('EXTNAME', 'LOGNU')
hdulist[-2].header.update('EXTNAME', 'WAVEMICR')
hdulist[-1].header.update('EXTNAME', 'EBLZ')
hdulist[0].header.update('NAME', 'Haardt & Madau 2012')
hdulist.writeto(filename)

#===========================================================================
#===========================================================================

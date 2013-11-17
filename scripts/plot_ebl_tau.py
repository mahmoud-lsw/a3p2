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
#import scipy.integrate
import matplotlib.pyplot as plt
#
import a3p2
#import a3p2.astro.cosmo
#import a3p2.tools.config
#import a3p2.ebl.model as eblmodel
#from a3p2.constants import *
import a3p2.ebl.measurements

#===========================================================================
# Functions

#===========================================================================
# Main
def plot_ebl_tau(filebase, redshifts=[.1, .5, 4.], yrange=[1E-2, 1E3], xrange=[1E-2, 1E2], plottitle='',
                 plottotal=False, loglevel='INFO') :

    # Time it!
    t_start = time.clock()

    # Configure logging
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    mr_filenames = []
    try :
        fb_a = eval(filebase)
        for fb in fb_a :
            mr_filenames += glob.glob(fb)
    except :
        mr_filenames = glob.glob(filebase)

    if len(mr_filenames) == 0 :
        logging.error('Could not find files matching {0}'.format(filebase))
        return

    plt.figure()
    if plottitle :
        plt.title(plottitle)

    #cmstep = 1. / float(len(ebl_models))
    cmfrac = 0.
    ls_array = ['--', '-.', ':']
    cm_array = [plt.cm.Greens, plt.cm.Blues, plt.cm.Purples]
    #cm_array = [plt.cm.Blues, plt.cm.Reds, plt.cm.Purples, plt.cm.Oranges]    
    lw = 2.

    plt_e_nsteps = 150
    e_a = 10. ** np.linspace(np.log10(xrange[0]), np.log10(xrange[-1]), plt_e_nsteps)
    tau_total = np.zeros((len(redshifts), plt_e_nsteps))

    plt.xlabel('Energy (TeV)')
    plt.ylabel(r'Optical depth $\tau$')
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    plt.xlim(xrange)
    plt.ylim(yrange)

    for fi, filename in enumerate(mr_filenames) :
        logging.info('Reading optical depth from file {0}'.format(filename))
        f = pyfits.open(filename)
        z = f['TAU'].data
        z[z == 0.] = 1E-42
        tau_lookup =  scipy.interpolate.RectBivariateSpline(
            x=f['TAUE'].data.field('LOGE'),
            y=f['TAUZ'].data.field('LOGZ'),
            z=np.log10(z),
            kx=1, ky=1
            )
        for zi, z in enumerate(redshifts) :
            label='{0} z={1:.1f}'.format(f[0].header['NAME'], z)
            tau = 10. ** tau_lookup(np.log10(e_a), np.log10(z))[:, 0]
            plt.plot(
                e_a, tau,
                color=cm_array[zi % len(cm_array)](.7), lw=lw,
                label=label,
                ls=ls_array[fi % len(ls_array)]
                )
            tau_total[zi] += tau
        f.close()

    if plottotal :
        for zi, z in enumerate(redshifts) :
            label='Total z={0:.1f}'.format(z)
            plt.plot(
                e_a, tau_total[zi],
                color=cm_array[zi % len(cm_array)](.7), #color='black',
                lw=lw, ls='-', #ls_array[zi % len(ls_array)],
                label=label
                )

    plt.legend(prop={'size': 'x-small'}, frameon=False, loc='upper left', ncol=2)

    metatxt =  r''
    plt.text(
        1.005, 1.0,
        datetime.datetime.today().strftime('%d-%h-%Y %H:%M') + metatxt,
        fontsize='xx-small', ha='left', va='top',
        transform=plt.gca().transAxes, rotation='vertical',
        weight='light', color='black'
        )
    plt.show()

#===========================================================================
# Main function
if __name__ == '__main__':
    # We should switch to argparse soon (python v2.7++)
    # http://docs.python.org/library/argparse.html#module-argparse
    import optparse
    parser = optparse.OptionParser(
        usage='%prog <file_base_name> [options]',
        description='Plots EBL optical depths'
    )
    parser.add_option(
        '-z','--redshifts',
        dest='redshifts',
        default='[.1, .5, 4.]',
        help='Redshifts to be plotted [default: %default].'
    )
    parser.add_option(
        '-y','--yrange',
        dest='yrange',
        default='[1E-2, 1E3]',
        help='Plot y-axis range in units optical depth [default: %default].'
    )
    parser.add_option(
        '-x','--xrange',
        dest='xrange',
        default='[1E-2, 1E2]',
        help='Plot x-axis range in units of TeV [default: %default].'
    )
    parser.add_option(
        '-t','--plot-title',
        dest='plottitle',
        default='',
        help='Title to be added the plot [default: %default].'
    )
    parser.add_option(
        '--plot-total',
        dest='plottotal',
        default=False,
        action='store_true',
        help='Plot total EBL [default: %default].'
    )
    parser.add_option(
        '-l','--log-level',
        dest='loglevel',
        default='INFO',
        help='Amount of logging e.g. DEBUG, INFO, WARNING, ERROR [default: %default].'
    )

    options, args = parser.parse_args()

    if len(args) is 1:
        plot_ebl_tau(
            filebase=args[0],
            redshifts=eval(options.redshifts),
            yrange=eval(options.yrange),
            xrange=eval(options.xrange),
            plottitle=options.plottitle,
            plottotal=options.plottotal
            )
    else :
        parser.print_help()

#===========================================================================
#===========================================================================

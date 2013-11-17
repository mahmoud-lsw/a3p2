#===========================================================================
# Imports

#import sys
import time
import logging
import os
#import gc
#import gzip
import datetime
#import json

import pyfits
import numpy as np
import scipy.interpolate
import scipy.integrate
import matplotlib.pyplot as plt

import a3p2
#import a3p2.astro.cosmo
#import a3p2.tools.config
import a3p2.ebl.model as eblmodel
from a3p2.constants import *
import a3p2.ebl.measurements

#===========================================================================
# Functions

#===========================================================================
# Main
def plot_ebl_model_sfr(config_file,
                       showHB2006=False, showBL2006=False, showT2007=False, showTS2009=False,
                       loglevel='INFO') :
    # Time it!
    t_start = time.clock()

    # Configure logging
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # Read in & calculate EBL models
    config_dict = eval(str(open(config_file).read()))
    ebl_models = []
    zmax_a = []
    for model_dict in config_dict['models'] :
        logging.info(50 * '=')
        if 'name' in model_dict.keys() :
            logging.info('Parsing model {0}'.format(model_dict['name']))
        ebl_models.append(eblmodel.EBLModel(eblmodel.ConfigEBLModel(model_dict)))
        logging.info('Model configuration:')
        ebl_models[-1].config.print_()
        zmax_a.append(ebl_models[-1].config.zintmax)

    zmax_a = np.array(zmax_a)

    plt.figure()
    # Chabrier IMF
    #sfr_scalef = 1. / 1.65 #.65
    # Salpeter IMF
    sfr_scalef = 1.
    zmax = 4.
    bpl = lambda p,x : np.where(x < p[0], p[1] * (x / p[0]) ** -p[2],  p[1] * (x / p[0]) ** -p[3])

    ls_array = ['-', '--', '-.', ':']
    cm_array = [plt.cm.Greens, plt.cm.Blues, plt.cm.Purples]
    lw = 1.5

    if showHB2006 :
        # Load and plot SFR data
        sfr_fit_3sigma = np.loadtxt(
            '/Users/mraue/Stuff/work/pop3paper/data/sfr-compilation/hopkins_beacom_2006_3sigma_sala.txt'
            )
        plt.fill(
            10. ** sfr_fit_3sigma[:, 0] - 1., sfr_fit_3sigma[:, 1] / .77 * sfr_scalef,
            color=plt.cm.Blues(.15),
            fill=True,# label=r'HB06, 3$\sigma$'
            )
        sfr_fit_1sigma = np.loadtxt(
            '/Users/mraue/Stuff/work/pop3paper/data/sfr-compilation/hopkins_beacom_2006_1sigma_sala.txt'
            )
        plt.fill(
            10. ** sfr_fit_1sigma[:, 0] - 1., sfr_fit_1sigma[:, 1] / .77 * sfr_scalef,
            color=plt.cm.Blues(.3),
            fill=True,# label=r'HB06, 1$\sigma$'
            )
    if showHB2006 and 1 :
        sfr_data = np.loadtxt('/Users/mraue/Stuff/work/pop3paper/data/sfr-compilation/sfh_compilation.dat')
        sfr_data_y = 10.**sfr_data[:,3]
        #plt.errorbar(x=sfr_data[:,0], y=sfr_data_y,
        #             xerr=[sfr_data[:,2],sfr_data[:,1]],
        #             yerr=[sfr_data_y - 10. ** (sfr_data[:,3] - sfr_data[:,5]), 10. ** (sfr_data[:,3] + sfr_data[:,4]) - sfr_data_y],
        #             fmt='.', zorder=0, color='.8')
        plt.errorbar(x=sfr_data[:,0], y=sfr_data_y * sfr_scalef,
                     xerr=[sfr_data[:,2],sfr_data[:,1]],
                     yerr=np.array([sfr_data_y - 10. ** (sfr_data[:,3] - sfr_data[:,5]),
                                    10. ** (sfr_data[:,3] + sfr_data[:,4]) - sfr_data_y]) * sfr_scalef,
                     fmt='.', color='.75', mec='.75', label='', capsize=0.)

    if 1 :
        # Plot Raue, Kneiske, Mazin 2009 limit
        limx = np.array([7., 15., 15., 7., 7.])
        limy = np.array([0.3, 0.75, 3., 1.1, 0.3])
        plt.fill(limx, limy, color=plt.cm.Purples(.3), fill=True, label='Raue et al. 2009')

    showB2012 = True 
    if showB2012 :
        # Plot Bouwens et al. 2012 SFRD measurements at higher redshifts
        z = np.array([3.8, 5., 5.9, 6.8, 8.])
        sfrd = np.array([-1.48, -1.78, -1.83, -2.024, -2.249])
        sfrd_err = np.array([0.05, 0.06, 0.08, 0.092, 0.1213])
        sfrd_cor = np.array([-1.05, -1.48, -1.65, -1.892, -2.128])
        plt.errorbar(
            x=z,
            y=10. ** sfrd * sfr_scalef,
            yerr=np.array([10. ** (sfrd + sfrd_err) - 10. ** sfrd, 10. ** sfrd - 10. ** (sfrd - sfrd_err)]) * sfr_scalef,
            label='B2012',
            fmt='o',
            color=(.2745, 0., 1.),
            mec=(.2745, 0., 1.),
            mfc='none',
            capsize=0.
            )
        plt.errorbar(
            x=z,
            y=10. ** sfrd_cor * sfr_scalef,
            yerr=np.array([10. ** (sfrd_cor + sfrd_err) - 10. ** sfrd_cor, 10. ** sfrd_cor - 10. ** (sfrd_cor - sfrd_err)]) * sfr_scalef,
            label='B2012_cor',
            fmt='o',
            c='r',
            mec='r',
            mfc='none',
            capsize=0.
            )

    if showT2007 :
        # Plot Tornatore et al. 2007 SFRD
        d = np.loadtxt(a3p2.datapath + '/sfr/tornatore_2007_fig1_pop2_10Mpc.txt')
        plt.plot(d[:,0], 10. ** d[:,1], '-', color='.6', lw=1.5, label='T2007 PopII')
        d = np.loadtxt(a3p2.datapath + '/sfr/tornatore_2007_fig1_pop3_10Mpc.txt')
        plt.plot(d[:,0], 10. ** d[:,1], '--', color='.6', lw=1.5, label='T2007 PopIII')

    if showTS2009 :
        # Plot Trenti & Stiavelli 2009 Fig. 6 SFRD
        d = np.loadtxt(a3p2.datapath + '/sfr/trenti_2009_fig6_pop2.txt')
        plt.plot(d[:,0], 10. ** d[:,1], '-', color='.5', lw=1.5, label='TS2009 PopII')
        d = np.loadtxt(a3p2.datapath + '/sfr/trenti_2009_fig6_pop3.txt')
        #plt.plot(d[:,0], 10. ** d[:,1], '-.', color='.5', lw=1.5, label='TS2009 PopIII')
        dintp = scipy.interpolate.UnivariateSpline(d[:,0], d[:,1], s=0, k=1)
        d = np.loadtxt(a3p2.datapath + '/sfr/trenti_2009_fig6_pop3h2.txt')
        #plt.plot(d[:,0], 10. ** d[:,1], ':', color='.5', lw=1.5, label='TS2009 PopIIIH2')
        plt.plot(d[:,0], 10. ** (d[:,1]) + 10. ** dintp(d[:,0]), '--', color='.5', lw=1.5, label='TS2009 PopIIITot')
    if showBL2006 :
        # Plot Bromm & Loeb 2006
        #d = np.loadtxt(a3p2.datapath + '/sfr/bromm_2006_fig1_strong_pop2.txt')
        #plt.plot(d[:,0], d[:,1], '-', color='.6', lw=1.5, label='BL2006 PopII')
        #d = np.loadtxt(a3p2.datapath + '/sfr/bromm_2006_fig1_strong_pop3.txt')
        #plt.plot(d[:,0], d[:,1], '--', color='.6', lw=1.5, label='BL2006 PopIII')
        d = np.loadtxt(a3p2.datapath + '/sfr/bromm_2006_fig1_weak_pop2.txt')
        plt.plot(d[:,0], d[:,1], '-', color='.4', lw=1.5, label='BL2006 PopII')
        d = np.loadtxt(a3p2.datapath + '/sfr/bromm_2006_fig1_weak_pop3.txt')
        plt.plot(d[:,0], d[:,1], '--', color='.4', lw=1.5, label='BL2006 PopIII')

    #z = np.linspace(0, zmax_a.max(), 100)
    z = np.linspace(1E-1, 35., 150)
    ls_a = ['--', '-.', ':']
    sfrd_sum = np.zeros(len(z))
    if len(ebl_models) :
        for mi, m in enumerate(ebl_models) :
            plt.plot(
                z, m.data['sfr'](z + 1.), lw=lw, color=cm_array[mi % len(cm_array)](.7),
                ls=ls_array[mi % len(ls_a)], label=m.config.name
                )
            sfrd_sum +=  m.data['sfr'](z + 1.)

            plt.plot(z, sfrd_sum, lw=1.5, color='black', ls='-', label='Total')

    #plt.xlim(p['zmin'], p['zmax'])
    plt.gca().set_yscale('log')
    #plt.gca().set_xscale('log')
    #plt.xlim(1E-1, 35.)
    #plt.ylim(1E-7, 5.)
    plt.xlim(1E-1, 16.)
    plt.ylim(1E-4, 10.)
    plt.xlabel('Redshift')
    plt.ylabel('SFRD (M$_\odot$ yr$^{-1}$ Mpc$^{-3}$)')

    lg = plt.legend(prop={'size': 10}, numpoints=1, loc='upper right', frameon=False, ncol=2)
    #lg.get_frame().set_linewidth(0)
    #plt.legend(prop={'size':11}, ncol=2, frameon=False)

    metatxt = 'sfr_scalef={0:.2f}'.format(sfr_scalef)
    plt.text(
        1.005, 1.0,
        datetime.datetime.today().strftime('%d-%h-%Y %H:%M, ') + metatxt,
        fontsize='xx-small', ha='left', va='top',
        transform=plt.gca().transAxes, rotation='vertical',
        weight='light', color='black'
        )
    if 'plottitle' in config_dict.keys() :
        plt.title(config_dict['plottitle'])
    plt.show()

#===========================================================================
# Main function
if __name__ == '__main__':
    # We should switch to argparse soon (python v2.7++)
    # http://docs.python.org/library/argparse.html#module-argparse
    import optparse
    parser = optparse.OptionParser(
        usage='%prog <config_file> [options]',
        description='Calculates EBL models.'
    )
    parser.add_option(
        '--show-HB2006',
        dest='showHB2006',
        default=False,
        action='store_true',
        help='Show Hopkins & Beacom 2006 SFRD compilation for comparison [default: %default].'
    )
    parser.add_option(
        '--show-BL2006',
        dest='showBL2006',
        default=False,
        action='store_true',
        help='Show Bromm & Loeb 2006 SFRD for comparison [default: %default].'
    )
    parser.add_option(
        '--show-T2007',
        dest='showT2007',
        default=False,
        action='store_true',
        help='Show Tornatore et al. 2007 SFRD for comparison [default: %default].'
    )
    parser.add_option(
        '--show-TS2009',
        dest='showTS2009',
        default=False,
        action='store_true',
        help='Show Trenti & Stavelli 2009 SFRD for comparison [default: %default].'
    )

    parser.add_option(
        '-l','--log-level',
        dest='loglevel',
        default='INFO',
        help='Amount of logging e.g. DEBUG, INFO, WARNING, ERROR [default: %default].'
    )

    options, args = parser.parse_args()

    if len(args) is 1:
        plot_ebl_model_sfr(
            config_file=args[0],
            showHB2006=options.showHB2006,
            showBL2006=options.showBL2006,
            showT2007=options.showT2007,
            showTS2009=options.showTS2009,
            loglevel=options.loglevel
            )
    else :
        parser.print_help()

#===========================================================================
#===========================================================================

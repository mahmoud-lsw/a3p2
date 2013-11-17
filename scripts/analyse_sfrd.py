#===========================================================================
# Imports

import time
import logging
import os
import datetime
import subprocess

import pyfits
import numpy as np
#import scipy.interpolate
import scipy.integrate
import scipy.special
import matplotlib.pyplot as plt

import ../a3p2
#import a3p2.astro.cosmo
from ../a3p2.tools.config import ConfigBase
import ../a3p2.ebl.model as eblmodel
from ../a3p2.constants import *
import ../a3p2.ebl.measurements

#===========================================================================
# Functions

#---------------------------------------------------------------------------
CONFIG_SCAN_SFRD_DEFAULT = {
    'beta_a' : [.3],
    'zpeakmin' : .2,
    'zpeakmax' : 2.4,
    'zpeakn' : 14,
    'lrho0min' : -1.2,
    'lrho0max' : -.5,
    'lrho0n' : 8,
    'rhoz0' : .012
    }

#---------------------------------------------------------------------------
class ConfigScanSFRD(ConfigBase) :
    def __init__(self, dict_=None, verbose=True) :
        ConfigBase.__init__(self, CONFIG_SCAN_SFRD_DEFAULT)
        if dict_ : ConfigBase.init_from_dict(self, dict_, verbose)

#===========================================================================
# Main
def analyse_sfrd(
        config_file,
        output_filename_base,
        do_graphical_output=True,
        loglevel='INFO'
        ) :

    # Time it!
    t_start = time.clock()

    # Configure logging
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # Read EBL models
    config_dict = eval(str(open(config_file).read()))
    ebl_models = []
    for model_dict in config_dict['models'] :
        logging.info(50 * '=')
        if 'name' in model_dict.keys() :
            logging.info('Parsing model {0}'.format(model_dict['name']))
        ebl_models.append(eblmodel.EBLModel(eblmodel.ConfigEBLModel(model_dict)))
        logging.info('Model configuration:')
        ebl_models[-1].config.print_()

    scan_sfrd_config = ConfigScanSFRD(config_dict['scansfrd'])
    scan_sfrd_config.print_()
    bpl_f = 'lambda p,x : np.where(x < p[0], p[1] * (x / p[0]) ** -p[2],  p[1] * (x / p[0]) ** -p[3])'
    bpl = lambda p,x : np.where(x < p[0], p[1] * (x / p[0]) ** -p[2],  p[1] * (x / p[0]) ** -p[3])
    # SFRD grid settings
    # SFR is a broken power law with first index being alpha and the second one beta
    beta_a = scan_sfrd_config.beta_a
    # SFR(z=0), Chab: 0.012, Salp: 0.02
    sfr_z0 = scan_sfrd_config.rhoz0
    #sfr_z0 = .02
    # GRID 1
    sfr_zpeak_a = np.linspace(scan_sfrd_config.zpeakmin, scan_sfrd_config.zpeakmax, scan_sfrd_config.zpeakn)
    sfr_l_a = np.linspace(scan_sfrd_config.lrho0min, scan_sfrd_config.lrho0max, scan_sfrd_config.lrho0n)

    neblmodels = len(sfr_zpeak_a) * len(sfr_l_a) * len(beta_a)
    logging.info('Total number of models to be evaluated: {0}'.format(neblmodels))

    results_ul_a = np.zeros([len(beta_a), len(sfr_zpeak_a), len(sfr_l_a)])
    #results_ll_a = np.zeros([len(beta_a), len(sfr_zpeak_a), len(sfr_l_a)])
    results_ll_a = np.ones([len(beta_a), len(sfr_zpeak_a), len(sfr_l_a)])

    if len(output_filename_base) > 5 and output_filename_base[-5:] == '.fits' :
        output_filename_base = output_filename_base[:-5]

    # Load EBL UL limit data
    ebl_limits = np.loadtxt(a3p2.datapath + '/ebl/limits/meyer_raue_et-al_2010_ul_final.dat')
    lebl_limits = np.log10(ebl_limits)

    plottitle = ''
    # Prepare EBL lower limits
    #ebl_ll_used = ['dole2006', 'elbaz2002', 'fazio2005', 'metcalfe2003', 'papovich2004',
    #               'bethermin2010ll', 'berta2010', 'madauPozzetti2000']
    # STD
    #plottitle = 'STD - '
    #ebl_ll_used = ['dole2006', 'elbaz2002', 'fazio2005', 'metcalfe2003', 'papovich2004',
    #               'bethermin2010ll', 'madauPozzetti2000']
    # EXTRA UV
    #plottitle = 'EXTRA UV - '
    #ebl_ll_used = ['dole2006', 'elbaz2002', 'fazio2005', 'metcalfe2003', 'papovich2004',
    #               'bethermin2010ll', 'madauPozzetti2000', 'xu2005']
    # NO FAZIO
    plottitle = 'NO FAZIO - '    
    ebl_ll_used = ['dole2006', 'elbaz2002', 'metcalfe2003', 'papovich2004',
                   'bethermin2010ll', 'madauPozzetti2000']
    #ebl_ll_used = ['fazio2005']
    #ebl_ll_used = ['madauPozzetti2000']
    #ebl_ll_used = ['berta2010']
    ebl_ll_col = a3p2.ebl.measurements.get_ebl_measurement_collection()
    ebl_ll_w, ebl_ll_val, ebl_ll_err = [], [], []
    for m in ebl_ll_col :
        if m['id'] in ebl_ll_used :
            logging.debug('Adding {0} to LL analysis'.format(m['id']))
            ebl_ll_w += m['data']['lambda']
            ebl_ll_val += m['data']['ebl']
            ebl_ll_err += m['data']['ebl_err_low']
    ebl_ll_w = np.array(ebl_ll_w)
    ebl_ll_val = np.array(ebl_ll_val)
    ebl_ll_err = np.array(ebl_ll_err)

    # Evaluate EBL vs limits
    nstep = 0
    for model in ebl_models :
        model.config.sfrf = bpl_f
        for i_beta, beta in enumerate(beta_a) :
            for i_sfr_zpeak, sfr_zpeak in enumerate(sfr_zpeak_a) :
                for i_sfr_l, sfr_l in enumerate(sfr_l_a) :
                    nstep += 1
                    logging.debug('::::::::::::::::: Step {0} of {1}'.format(nstep, neblmodels))
                    logging.debug('beta = {0}, zpeak = {1}, lrho0 = {2}'.format(beta, sfr_zpeak, sfr_l))
                    alpha = (sfr_l - np.log10(sfr_z0)) / (np.log10(sfr_zpeak + 1.) - np.log10(0. + 1.))
                    sfrp0 = [sfr_zpeak + 1., 10. ** sfr_l, -alpha, beta]
                    model.config.sfrp0 = sfrp0
                    sfr = lambda x: bpl(sfrp0, x) # [M_solar Mpc^-3]
                    logging.debug('SFR : {0}'.format(sfrp0))
                    filename = output_filename_base + '_' + model.config.name.replace(' ', '_')
                    filename += '_n{0}.model.fits.gz'.format(nstep)
                    logging.debug('Reading EBL model from file {0}'.format(filename))
                    f = pyfits.open(filename)
                    h = f[0].header
                    # Check if data file matches config file
                    san_d = {'alpha' : alpha, 'beta' : beta, 'lrho0': sfr_l, 'zpeak': sfr_zpeak}
                    for key, val in san_d.iteritems() :
                        if np.abs(val - h[key]) > 1E-8 :
                            logging.error('{0} values do not match: {1} vs {2}'.format(key, val, h[key]))
                            return
                    d_z = f['EBLZ'].data.field('Z')
                    d_l = f['WAVEMICR'].data.field('WAVE')
                    d_ebl = f['EBL'].data
                    ebl_intp_z0 = scipy.interpolate.UnivariateSpline(
                        x = np.log10(d_l[::-1]),
                        y = np.log10(d_ebl[:, 0][::-1]),
                        s=0, k=1
                        )
                    # Handle upper limits
                    ul_frac = 10. ** (ebl_intp_z0(lebl_limits[:,0]) - lebl_limits[:,1])
                    results_ul_a[i_beta, i_sfr_zpeak, i_sfr_l] = ul_frac.max()
                    logging.debug('ul_frac.max() = {0}'.format(ul_frac.max()))
                    # Handle lower limits
                    ll_delta = 10. ** ebl_intp_z0(np.log10(ebl_ll_w)) - ebl_ll_val
                    m = ll_delta < 0.
                    ll_chi2 = np.sum((ll_delta[m] / ebl_ll_err[m]) ** 2.)
                    if ll_chi2 > 0. :
                        #results_ll_a[i_beta, i_sfr_zpeak, i_sfr_l] = ll_chi2 / np.sum(m)
                        results_ll_a[i_beta, i_sfr_zpeak, i_sfr_l] = np.log10(1. - scipy.special.gammainc(.5 * np.sum(m), .5 * ll_chi2))
                    logging.debug(
                        'll_chi / d.o.f. = {0} / {1} = {2}'.format(ll_chi2, np.sum(m), ll_chi2 / np.sum(m))
                        )
                    f.close()

    # Time it!
    logging.info('Excecution took {0:.2f} s'.format(time.clock() - t_start))
    
    # Graphical output
    if do_graphical_output :
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

        showHB2006 = True
        
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
        if showHB2006 and 0 :
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
                         fmt='.', color='.75', mec='.75')

        if 0 :
            # Plot Raue, Kneiske, Mazin 2009 limit
            limx = np.array([7., 15., 15., 7., 7.])
            limy = np.array([0.3, 0.75, 3., 1.1, 0.3])
            plt.fill(limx, limy, color=plt.cm.Purples(.3), fill=True, label='Raue et al. 2009')

        for i, beta in enumerate([.3]) :#(beta_a) :
            CS = plt.contour(
                sfr_zpeak_a, 10. ** sfr_l_a, np.transpose(results_ul_a[i,:,:]), [1., 1.2, 1.5],
                colors=plt.cm.Reds(.6), label='beta={0:.1f}'.format(beta),
                linestyles=ls_array[i % len(ls_array)], linewidths=lw
                )
            plt.clabel(CS, inline=1, fontsize=8, fmt='%1.1f')
            #CS = plt.contour(
            #    sfr_zpeak_a, 10. ** sfr_l_a, np.transpose(results_ll_a[i,:,:]), [2., 5., 10.],
            #    colors=plt.cm.Blues(.6), label='beta={0:.1f}'.format(beta),
            #    linestyles=ls_array[i % len(ls_array)], linewidths=lw
            #    )
            #plt.clabel(CS, inline=1, fontsize=8, fmt='%1.0f')
            CS = plt.contour(
                sfr_zpeak_a, 10. ** sfr_l_a, np.transpose(results_ll_a[i,:,:]), [-2., -4., -6.],
                colors=plt.cm.Blues(.6), label='beta={0:.1f}'.format(beta),
                linestyles=ls_array[i % len(ls_array)], linewidths=lw
                )
            plt.clabel(CS, inline=1, fontsize=8, fmt='%1.0f')

        #plt.xlim(p['zmin'], p['zmax'])
        plt.gca().set_yscale('log')
        #plt.gca().set_xscale('log')
        plt.xlim(0., 3.)
        plt.ylim(1E-2, .5)
        plt.xlabel('Redshift')
        plt.ylabel('SFRD (M$_\odot$ yr$^{-1}$ Mpc$^{-3}$)')

        lg = plt.legend(prop={'size': 10}, numpoints=1, loc='upper right', frameon=False, ncol=2)

        metatxt = r''#', {0}'.format(str(ebl_ll_used))
        plt.text(
            1.005, 1.0,
            datetime.datetime.today().strftime('%d-%h-%Y %H:%M') + metatxt,
            fontsize=8, ha='left', va='top',
            transform=plt.gca().transAxes, rotation='vertical',
            weight='light', color='black'
            )
        if 'plottitle' in config_dict.keys() :
            plottitle += config_dict['plottitle']
        plt.title(plottitle)
        plt.show()

#===========================================================================
# Main function
if __name__ == '__main__':
    # We should switch to argparse soon (python v2.7++)
    # http://docs.python.org/library/argparse.html#module-argparse
    import optparse
    parser = optparse.OptionParser(
        usage='%prog <config_file> <output_filename_base> [options]',
        description='Calculates EBL models for a set of SFRD parameters.'
    )
    parser.add_option(
        '--no-graphical-output',
        dest='graphical_output',
        action='store_false',
        default=True,
        help='Switch off graphical output.'
    )
    parser.add_option(
        '-l','--log-level',
        dest='loglevel',
        default='INFO',
        help='Amount of logging e.g. DEBUG, INFO, WARNING, ERROR [default: %default].'
    )

    options, args = parser.parse_args()

    if len(args) is 2:
        analyse_sfrd(
            config_file=args[0],
            output_filename_base=args[1],
            do_graphical_output=options.graphical_output,
            loglevel=options.loglevel
            )
    else :
        parser.print_help()

#===========================================================================
#===========================================================================

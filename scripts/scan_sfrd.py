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
import matplotlib.pyplot as plt

import a3p2
#import a3p2.astro.cosmo
from a3p2.tools.config import ConfigBase
import a3p2.ebl.model as eblmodel
from a3p2.constants import *
import a3p2.ebl.measurements

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
def scan_sfrd(
        config_file,
        output_filename_base=None,
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
    logging.info('Total number of models to be calculated: {0}'.format(neblmodels))
    logging.info('Estimated calculation time : {0:.1f} min'.format(neblmodels * .45))

    if output_filename_basena and len(output_filename_base) > 5 and output_filename_base[-5:] == '.fits' :
        output_filename_base = output_filename_base[:-5]

    # Calculate EBL
    nstep, fsize_total = 0, 0.
    for model in ebl_models :
        model.config.sfrf = bpl_f
        for i_beta, beta in enumerate(beta_a) :
            for i_sfr_zpeak, sfr_zpeak in enumerate(sfr_zpeak_a) :
                for i_sfr_l, sfr_l in enumerate(sfr_l_a) :
                    nstep += 1
                    logging.info('::::::::::::::::: Step {0} of {1}'.format(nstep, neblmodels))
                    logging.info('beta = {0}, zpeak = {1}, lrho0 = {2}'.format(beta, sfr_zpeak, sfr_l))
                    alpha = (sfr_l - np.log10(sfr_z0)) / (np.log10(sfr_zpeak + 1.) - np.log10(0. + 1.))
                    sfrp0 = [sfr_zpeak + 1., 10. ** sfr_l, -alpha, beta]
                    model.config.sfrp0 = sfrp0
                    sfr = lambda x: bpl(sfrp0, x) # [M_solar Mpc^-3]
                    logging.debug('SFR : {0}'.format(sfrp0))
                    model.data['sfr'] = sfr
                    model.calculate(fesc=model.config.fesc)
                    # Write EBL model data to file
                    if output_filename_base :
                        filename = output_filename_base + '_' + model.config.name.replace(' ', '_')
                        filename += '_n{0}.model.fits'.format(nstep)
                        logging.info('Writing EBL model to file {0}'.format(filename))
                        model.write_fits(filename)
                        f = pyfits.open(filename, 'update')
                        scan_sfrd_config.write_to_fits_header(f[0].header)
                        f[0].header.update('beta', beta)
                        f[0].header.update('alpha', alpha)
                        f[0].header.update('zpeak', sfr_zpeak)
                        f[0].header.update('lrho0', sfr_l)
                        f[0].header.update('sfrp0', str(sfrp0))
                        f[0].header.update('ibeta', i_beta)
                        f[0].header.update('izpeak', i_sfr_zpeak)
                        f[0].header.update('irho0', i_sfr_l)
                        f.close()
                        logging.info('Compressing file ...')
                        subprocess.call(['gzip', filename])
                        fsize = os.path.getsize(filename + '.gz') / 1024.
                        logging.info('Written {0:.0f} KB.'.format(fsize))
                        fsize_total += fsize

    # Time it!
    logging.info('Excecution took {0:.2f} s'.format(time.clock() - t_start))
    logging.info('Written {0:.2f} MB total.'.format(fsize_total / 1024.))

    # Graphical output
    if do_graphical_output :
        plt.figure()
        # Plot EBL measurements and limits
        a3p2.ebl.measurements.plot_ebl_measurement_collection(nolabel=True)
        plt_ebl_nsteps = 150
        ebl_lambda = 10. ** np.linspace(np.log10(.1), np.log10(1E3), plt_ebl_nsteps)
        ebl_em_total = np.zeros(plt_ebl_nsteps)
        
        plt.xlabel(r'Wavelength ($\mu$m)')
        plt.ylabel(r'EBL SED (nW / m$^2$ sr)')
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.xlim([.1, 1E3])
        plt.ylim([1E-6, 1E2])

        #plt.legend(prop={'size':11}, ncol=2, frameon=False)

        metatxt =  r''
        plt.text(
            1.01, .0,
            datetime.datetime.today().strftime('%d-%h-%Y %H:%M') + metatxt,
            fontsize=8, ha='left', va='bottom',
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
        description='Calculates EBL models for a set of SFRD parameters.'
    )
    parser.add_option(
        '-o','--output-filename-base',
        dest='output_filename_base',
        type='string',
        default=None,
        help='Output filename base.'
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

    if len(args) is 1:
        scan_sfrd(
            config_file=args[0],
            output_filename_base=options.output_filename_base,
            do_graphical_output=options.graphical_output,
            loglevel=options.loglevel
            )
    else :
        parser.print_help()

#===========================================================================
#===========================================================================

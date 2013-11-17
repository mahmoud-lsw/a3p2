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

#---------------------------------------------------------------------------
CONFIG_EBL_ATTENUATION_DEFAULT = {
    'name' : 'default',
    'zmin' : .001,
    'zmax' : 8.,
    'lznsteps' : 61,
    'emin' : 0.01,
    'emax' : 100.,
    'lensteps' : 61,
    'cosmo' : [.7, .3, .7]
    }

#---------------------------------------------------------------------------
class ConfigEBLAttenuation(ConfigBase) :
    def __init__(self, dict_=None, verbose=True) :
        ConfigBase.__init__(self, CONFIG_EBL_ATTENUATION_DEFAULT)
        if dict_ : ConfigBase.init_from_dict(self, dict_, verbose)

#===========================================================================
# Main
def calculate_ebl_attenuation(model_file,
                              config_file,
                              add_tau=False,
                              output_filename=None,
                              do_graphical_output=True,
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

    # Read configuration
    config = ConfigEBLAttenuation()
    if config_file :
        logging.debug('Opening configuration file {0}'.format(config_file))
        config = ConfigEBLAttenuation(eval(str(open(config_file).read())))
    else :
        logging.debug('Using default configuration')
    config.print_()

    # Read in EBL model
    logging.debug('Opening EBL model file {0}'.format(model_file))
    ebl = EBLModelMR2012(model_file)

    le_a = np.linspace(np.log10(config.emin), np.log10(config.emax), config.lensteps)
    lz_a = np.linspace(np.log10(config.zmin), np.log10(config.zmax), config.lznsteps)
    tau = np.zeros([config.lensteps, config.lznsteps])
    
    logging.info('Calculating optical depth ...')

    for ei, e in enumerate(le_a) :
        logging.debug('Step {0} of {1}'.format(ei + 1, len(le_a)))
        for zi, z in enumerate(lz_a) :
            tau[ei, zi] = calc_ebl_attenuation2(
                10. ** z, 10. ** e, ebl,
                h=config.cosmo[0], W_m=config.cosmo[1], W_l=config.cosmo[2]
                )

    def tau_to_hdulist(tau, e_a, z_a, hdul=None) :
        if hdul == None :
            hdul = pyfits.HDUList([pyfits.PrimaryHDU()])
        hdul.append(pyfits.ImageHDU(data=tau, name='TAU'))
        hdul.append(pyfits.new_table(pyfits.ColDefs([
            pyfits.Column(name='LOGE',format='E', array=le_a, unit='log10(TeV)')
            ])),)
        hdul.append(pyfits.new_table(pyfits.ColDefs([
            pyfits.Column(name='LOGZ',format='E', array=lz_a)
            ])),)
        hdul[-2].header.update('EXTNAME', 'TAUE')
        hdul[-1].header.update('EXTNAME', 'TAUZ')
        return hdul

    # Write EBL attenuation to file
    if output_filename :
        logging.info('Writing tau data to file {0}'.format(output_filename))
        f = tau_to_hdulist(tau, le_a, lz_a)
        config.write_to_fits_header(f['TAU'].header)
        #f[0].header.update('name', ebl.name)
        # Copy EBL model header
        fmodel = pyfits.open(model_file)
        for card in fmodel[0].header.ascardlist() :
            f[0].header.update(card.key, card.value, card.comment)
        fmodel.close()
        f.writeto(output_filename)
        logging.info('Written {0:.0f} KB.'.format(os.path.getsize(output_filename)/1024.))

    # Add EBL attenuation to EBL model file
    if add_tau :
        logging.info('Adding tau data to file {0}'.format(model_file))
        f = pyfits.open(model_file, mode='update')
        tau_to_hdulist(tau, le_a, lz_a, f)
        config.write_to_fits_header(f['TAU'].header)
        f.flush()
        logging.info('Done.')
    
    # Time it!
    logging.info('Excecution took {0:.2f} s'.format(time.clock() - t_start))

    # Graphical output
    if do_graphical_output :
        plt.figure()
        for zi in range(len(lz_a)) :
            plt.loglog(10. ** le_a, tau[:, zi], label='z={0:.3f}'.format(10. ** lz_a[zi]))
        plt.legend(loc='lower right', prop={'size' : 10}, frameon=False, ncol=3)
        plt.xlim([config.emin, config.emax])
        plt.ylim([1E-5, tau.max()])
        plt.show()
    
#===========================================================================
# Main function
if __name__ == '__main__':
    # We should switch to argparse soon (python v2.7++)
    # http://docs.python.org/library/argparse.html#module-argparse
    import optparse
    parser = optparse.OptionParser(
        usage='%prog <ebl_model_file> [options]',
        description='Calculates EBL model attenuation.'
    )
    parser.add_option(
        '-c','--config-filename',
        dest='config_filename',
        type='string',
        default=None,
        help='Output filename.'
    )
    parser.add_option(
        '-o','--output-filename',
        dest='output_filename',
        type='string',
        default=None,
        help='Output filename.'
    )
    parser.add_option(
        '-a', '--add-tau-to-model-file',
        dest='add_tau',
        action='store_true',
        default=False,
        help='Add optical depth to EBL model file.'
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
        calculate_ebl_attenuation(
            model_file=args[0],
            config_file=options.config_filename,
            output_filename=options.output_filename,
            add_tau=options.add_tau,
            do_graphical_output=options.graphical_output,
            loglevel=options.loglevel
            )
    else :
        parser.print_help()

#===========================================================================
#===========================================================================

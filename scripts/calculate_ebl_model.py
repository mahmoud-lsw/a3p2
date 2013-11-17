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
#import scipy.interpolate
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
def calculate_ebl_model(config_file,
                        output_filename_base=None,
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

    # Read in & calculate EBL models
    config_dict = eval(str(open(config_file).read()))
    ebl_models = []
    for model_dict in config_dict['models'] :
        logging.info(50 * '=')
        if 'name' in model_dict.keys() :
            logging.info('Parsing model {0}'.format(model_dict['name']))
        ebl_models.append(eblmodel.EBLModel(eblmodel.ConfigEBLModel(model_dict)))
        logging.info('Model configuration:')
        ebl_models[-1].config.print_()
        logging.info('Calculating model ...')
        ebl_models[-1].calculate(fesc=ebl_models[-1].config.fesc)
        logging.info('Calculation finished.')

    # Write EBL model data to file
    if output_filename_base :
        if len(output_filename_base) > 5 and output_filename_base[-5:] == '.fits' :
            output_filename_base = output_filename_base[:-5]
        for model in ebl_models :
            filename = output_filename_base + '_' + model.config.name.replace(' ', '_') + '.model.fits'
            logging.info('Writing EBL model {0} to file {1}'.format(model.config.name, filename))
            model.write_fits(filename)
            logging.info('Written {0:.0f} KB.'.format(os.path.getsize(filename)/1024.))

    # Time it!
    logging.info('Excecution took {0:.2f} s'.format(time.clock() - t_start))

    # Graphical output
    if do_graphical_output :
        plt.figure()
        # Plot EBL measurements and limits
        a3p2.ebl.measurements.plot_ebl_measurement_collection(nolabel=True)
        plt_ebl_nsteps = 150
        ebl_lambda = 10. ** np.linspace(np.log10(.1), np.log10(1E3), plt_ebl_nsteps)
        ebl_em_total = np.zeros(plt_ebl_nsteps)

        cmstep = 1. / float(len(ebl_models))
        cmfrac = 0.
        ls_array = ['-.', '--', '-', ':']
        for lsx, model in enumerate(ebl_models) :
            ebl_em = 10. ** model.data['ebl_intp'].ev(
                np.log10(si.c / ebl_lambda / 1E-6),
                np.zeros(len(ebl_lambda))
                )
            if lsx > len(ls_array) - 1 :
                lsx = lsx % len(ls_array)
            plt.loglog(
                ebl_lambda, ebl_em,
                ls=ls_array[lsx], color='black', label=r'{0}'.format(model.config.name), lw=2.
                )
            ebl_em_total += ebl_em
            cmfrac += cmstep
        plt.loglog(
            ebl_lambda, ebl_em_total,
            ls='-', color='black', label=r'Total', lw=2.
            )

        # DEBUG
        #linti = (np.abs(ebl_lambda - 5.5)).argmin()
        #lintmaxi = (np.abs(ebl_lambda - 1E3)).argmin()
        #print scipy.integrate.simps(y=ebl_em_total[linti:lintmaxi] * ebl_lambda[linti:lintmaxi] + np.log(10.),
        #                            x=np.log10(ebl_lambda[linti:lintmaxi]))
        # DEBUG END

        plt.xlabel(r'Wavelength ($\mu$m)')
        plt.ylabel(r'EBL SED (nW / m$^2$ sr)')
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.xlim([.1, 1E3])
        plt.ylim([1E-6, 1E2])

        pathtofile =  a3p2.datapath + '/ebl/models/test/'
        #kd2010datafiles = [
        #    'BC_salp02_MRF_Opt_1.20_0.07_3.50_-1.20_0.05_IR_1.20_0.06_4.50_-1.00_0.40.dat',
        #    'BC_salp02_MRF_IR_1.20_0.07_3.50_-1.20_0.05_IR_1.20_0.06_4.50_-1.00_0.40.dat',
        #    'BC_salp02_MRF_Tot_1.20_0.07_3.50_-1.20_0.05_IR_1.20_0.06_4.50_-1.00_0.40.dat'
        #    ]
        kd2010datafiles = [
            'BC_salp02_MRF_Opt_1.20_0.07_3.50_-1.20_0.05_IR_1.20_0.06_4.50_-1.00_0.80.dat',
            'BC_salp02_MRF_IR_1.20_0.07_3.50_-1.20_0.05_IR_1.20_0.06_4.50_-1.00_0.80.dat',
            'BC_salp02_MRF_Tot_1.20_0.07_3.50_-1.20_0.05_IR_1.20_0.06_4.50_-1.00_0.80.dat'
            ]
        for lsx, modelfile in enumerate(kd2010datafiles) :
            ebl_model_data = np.loadtxt(pathtofile + modelfile, skiprows=1)
            modelid = modelfile.split('_')[3]
            plt.plot(ebl_model_data[:,0], ebl_model_data[:,1], color=plt.cm.Blues(.7), lw=2.,
                     ls=ls_array[lsx], label='KD2010 ' + modelid)
            # DEBUG
            #linti = (np.abs(ebl_model_data[:,0] - 5.5)).argmin()
            #lintmaxi = (np.abs(ebl_model_data[:,0] - 1E3)).argmin()
            #print scipy.integrate.simps(
            #    y=ebl_model_data[linti:lintmaxi, 1] * ebl_model_data[linti:lintmaxi,0] + np.log(10.),
            #    x=np.log10(ebl_model_data[linti:lintmaxi, 0])
            #    )
            # DEBUG END

        ebl_model_file = a3p2.datapath + '/ebl/models/kneiske_dole_2010_MRF_LL.dat'
        ebl_model_data = np.loadtxt(ebl_model_file, skiprows=1)
        plt.plot(ebl_model_data[:,0], ebl_model_data[:,1], color=plt.cm.Greens(.7), lw=2.,
                 ls='-', label='KD2010 - paper')

        plt.legend(prop={'size':11}, ncol=2, frameon=False)

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
        description='Calculates EBL models.'
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
        calculate_ebl_model(
            config_file=args[0],
            output_filename_base=options.output_filename_base,
            do_graphical_output=options.graphical_output,
            loglevel=options.loglevel
            )
    else :
        parser.print_help()

#===========================================================================
#===========================================================================

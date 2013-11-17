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
def plot_ebl_sed(filebase, redshifts, yrange=[1E-6, 1E2], xrange=[.1, 1E3], plottitle='',
                 plottotal=False, plotlegend=True, showkd2010=False, showhm2012=False,
                 showid=False, loglevel='INFO') :

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
    # Plot EBL measurements and limits
    a3p2.ebl.measurements.plot_ebl_measurement_collection(nolabel=True)


    #cmstep = 1. / float(len(ebl_models))
    cmfrac = 0.
    ls_array = ['--', '-.', ':']
    cm_array = [plt.cm.Greens, plt.cm.Blues, plt.cm.Purples]
    #cm_array = [plt.cm.Blues, plt.cm.Reds, plt.cm.Purples, plt.cm.Oranges]    
    lw = 2.

    redshifts = eval(redshifts)

    plt_ebl_nsteps = 150
    ebl_lambda = 10. ** np.linspace(np.log10(xrange[0]), np.log10(xrange[-1]), plt_ebl_nsteps)
    ebl_em_total = np.zeros((len(redshifts), plt_ebl_nsteps))

    #for lsx, model in enumerate(ebl_models) :
    #    ebl_em = 10. ** model.data['ebl_intp'].ev(
    #        np.log10(si.c / ebl_lambda / 1E-6),
    #        np.zeros(len(ebl_lambda))
    #        )
    #    if lsx > len(ls_array) - 1 :
    #        lsx = lsx % len(ls_array)
    #    plt.loglog(
    #        ebl_lambda, ebl_em,
    #        ls=ls_array[lsx], color='black', label=r'{0}'.format(model.config.name), lw=2.
    #        )
    #    ebl_em_total += ebl_em
    #    cmfrac += cmstep
    #plt.loglog(
    #    ebl_lambda, ebl_em_total,
    #    ls='-', color='black', label=r'Total', lw=2.
    #    )

    plt.xlabel(r'Wavelength ($\mu$m)')
    plt.ylabel(r'EBL SED (nW / m$^2$ sr)')
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    plt.xlim(xrange)
    plt.ylim(yrange)

    if  0 :
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
            #plt.plot(ebl_model_data[:,0], ebl_model_data[:,1], color=plt.cm.Blues(.7), lw=1.,
            #         ls=ls_array[lsx], label='KD2010 ' + modelid)


    if showkd2010 :
        ebl_model_file = a3p2.datapath + '/ebl/models/kneiske_dole_2010_MRF_LL.dat'
        ebl_model_data = np.loadtxt(ebl_model_file)
        for zi, z in enumerate(redshifts) :
            idx = (np.abs(ebl_model_data[0, 1:] - z)).argmin()
            if idx > -1 and idx < len(ebl_model_data[1:,0]) - 1:
                label = 'KD2010 z={0:.1f}'.format(z)
                if showid :
                    label = 'KD2010 z={0:.1f} id={1}'.format(z, idx)
                plt.plot(ebl_model_data[1:,0], ebl_model_data[1:, idx + 1], color='.5', lw=lw,
                         ls=ls_array[zi % len(ls_array)], label=label)
            else :
                logging.warning('Could not find matching KD2010 model for z={0}, idx={1}'.format(z, idx))

    if showhm2012 :
        ebl_model_file = a3p2.datapath + '/ebl/models/haardt_madau_2012_UVB.out'
        ebl_model_data = np.loadtxt(ebl_model_file)
        for zi, z in enumerate(redshifts) :
            idx = (np.abs(ebl_model_data[0, 1:] - z)).argmin()
            if idx > -1 and idx < len(ebl_model_data[1:,0]) - 1:
                label = 'HM2012 z={0:.1f}'.format(z)
                if showid :
                    label = 'HM2012 z={0:.1f} id={1}'.format(z, idx)
                l = ebl_model_data[1:,0] * 1e-4
                plt.plot(l,
                         ebl_model_data[1:, idx + 1] * 3000. / l * 1E17 * (1. + z) ** -3.,# ergs/s/cm^2/Hz/sr to nW m^-2 sr^-1
                         color='.6', lw=lw,
                         ls=ls_array[zi % len(ls_array)], label=label)
            else :
                logging.warning('Could not find matching HM2012 model for z={0}, idx={1}'.format(z, idx))

    for fi, filename in enumerate(mr_filenames) :
        f = pyfits.open(filename)
        mr_z = f['EBLZ'].data.field('Z')
        mr_l = f['WAVEMICR'].data.field('WAVE')
        mr_ebl = f['EBL'].data
        for zi, z in enumerate(redshifts) :
            idx = (np.abs(mr_z - z)).argmin()
            label='{0} z={1:.1f}'.format(f[0].header['NAME'], z)
            if showid :
                label='{0} z={1:.1f} id={2}'.format(f[0].header['NAME'], z, idx)
            plt.plot(
                mr_l[::-1], mr_ebl[:, idx][::-1], color=cm_array[zi % len(cm_array)](.7), lw=lw,
                label=label,
                ls=ls_array[fi % len(ls_array)]
                )
            ebl_em_total[zi] += 10. ** scipy.interpolate.UnivariateSpline(
                x=np.log10(mr_l[::-1]),
                y=np.log10(mr_ebl[:, idx][::-1]),
                s=0, k=1,
                )(np.log10(ebl_lambda))
        f.close()

    if plottotal :
        for zi, z in enumerate(redshifts) :
            label='Total z={0:.1f}'.format(z)
            plt.plot(
                ebl_lambda, ebl_em_total[zi],
                color=cm_array[zi % len(cm_array)](.7), #color='black',
                lw=lw, ls='-', #ls=ls_array[zi % len(ls_array)],
                label=label
                )

    if plotlegend :
        plt.legend(prop={'size': 'x-small'}, frameon=False, loc='lower right', ncol=2)

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
        description='Plots EBL SEDs.'
    )
    parser.add_option(
        '-z','--redshifts',
        dest='redshifts',
        default='[0., 2., 4.]',
        help='Redshifts to be plotted [default: %default].'
    )
    parser.add_option(
        '-y','--yrange',
        dest='yrange',
        default='[1E-6, 1E2]',
        help='Plot y-axis range in units of nW sr^-1 m^-2 [default: %default].'
    )
    parser.add_option(
        '-x','--xrange',
        dest='xrange',
        default='[.1, 1E3]',
        help='Plot x-axis range in units of micrometer [default: %default].'
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
        '--no-legend',
        dest='plotlegend',
        default=True,
        action='store_false',
        help='Do not add a legend to plot.'
    )
    parser.add_option(
        '--show-id',
        dest='showid',
        default=False,
        action='store_true',
        help='Show column number [default: %default].'
    )
    parser.add_option(
        '--show-kd2010',
        dest='showkd2010',
        default=False,
        action='store_true',
        help='Show Kneiske & Dole 2010 for comparison [default: %default].'
    )
    parser.add_option(
        '--show-hm2012',
        dest='showhm2012',
        default=False,
        action='store_true',
        help='Show Haart & Madau 2012 for comparison [default: %default].'
    )
    parser.add_option(
        '-l','--log-level',
        dest='loglevel',
        default='INFO',
        help='Amount of logging e.g. DEBUG, INFO, WARNING, ERROR [default: %default].'
    )

    options, args = parser.parse_args()

    if len(args) is 1:
        plot_ebl_sed(
            filebase=args[0],
            redshifts=options.redshifts,
            yrange=eval(options.yrange),
            xrange=eval(options.xrange),
            plottitle=options.plottitle,
            plottotal=options.plottotal,
            plotlegend=options.plotlegend,
            showkd2010=options.showkd2010,
            showhm2012=options.showhm2012,
            showid=options.showid,
            loglevel=options.loglevel
            )
    else :
        parser.print_help()

#===========================================================================
#===========================================================================

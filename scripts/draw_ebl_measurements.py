from __future__ import  absolute_import
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import a3p2
import a3p2.ebl.measurements as ebl

plt.figure(1, figsize=(9,6))

yearcut=None#2004
ebl.plot_ebl_measurement_collection(yearmax=yearcut, nolabel=True)
ebl.plot_ebl_measurement_collection(cm=plt.cm.Dark2, yearmin=yearcut)

plt.axis([0.1, .7E3, 0.2, .7E3])

plt.subplots_adjust(left=0.125, right=.725, top=.95, bottom=.13)

# Legend handler for errobar with marker only
my_handler = matplotlib.legend_handler.HandlerErrorbar(xerr_size=.001, yerr_size=None, marker_pad=.0)

plt.legend(numpoints=1,prop={'size':6}, bbox_to_anchor=(1.05, 1), loc=2, ncol=1,
           borderaxespad=0., frameon=False, handler_map={matplotlib.axes.ErrorbarContainer:my_handler})

plt.xlabel(r'Wavelength ($\mu$m)')
plt.ylabel(r'EBL SED (nW / m$^2$ sr)')

plt.gca().set_yscale('log')
plt.gca().set_xscale('log')

plt.show()

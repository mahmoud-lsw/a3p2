#===========================================================================
# Imports

import logging

import numpy as np
import scipy.optimize
import scipy.special

#===========================================================================
# Functions & classes

#---------------------------------------------------------------------------
class ChisquareFitter :
    """
    Convenience class to perform Chi^2 fits.

    Attributes
    ----------

    fitfunc : function
        Fitfunction
    results : array
        Fit results from scipy.optimize.leastsq()
    chi_arr : float array
        Array with the final chi values
    chi2 : float
        Summed Chi^2
    dof : float
        Degrees of freedom
    prob : float
        Probability of the fit

    Parameters
    ----------

    fitfunc : function
        Fit function.
    """

    def __init__(self, fitfunc) :
        self.fitfunc = fitfunc
        self.results = None


    def fit_data(self, p0, x, y, y_err) :
        """
        Perform actual fit.

        Parameters
        ----------

        p0 : float array
            Start parameters
        x, y, y_err : float arrays
            Data to be fitted.
        """
        self.results = scipy.optimize.leastsq(self.chi_func, p0, args=(x, y, y_err), full_output=True)
        if self.results[4] :
            self.chi_arr = self.chi_func(self.results[0], x, y, y_err)
            self.chi2 = np.sum(np.power(self.chi_arr, 2.))
            self.dof = len(x) - len(p0)
            #self.prob = scipy.special.gammainc(.5 * self.dof, .5 * self.chi2) / scipy.special.gamma(.5 * self.dof)
            self.prob = 1. - scipy.special.gammainc(.5 * self.dof, .5 * self.chi2)
        return self.results[4]

    def chi_func(self, p, x, y, err):
        """Returns Chi"""
        return (self.fitfunc(p, x) - y) / err # Distance to the target function

    def print_results(self) :
        """Prints out results to the command line using the logging module."""
        if self.results == None :
            logging.warning('No fit results to report since no fit has been performed yet')
            return
        if self.results[4] < 5 :
            logging.info('Fit was successful!')
        else :
            logging.warning('Fitting failed!')
            logging.warning('Message: {0}'.format(self.results[3]))
        logging.info('Chi^2  : {0:f}'.format(self.chi2))
        logging.info('d.o.f. : {0:d}'.format(self.dof))
        logging.info('Prob.  : {0:.4e}'.format(self.prob))
        for i, v in enumerate(self.results[0]) :
            if self.results[1] != None :
                logging.info('P{0}     : {1:.4e} +/- {2:.4e}'.format(i, v,
                                                                     np.sqrt(self.results[1][i][i])))
            else :
                logging.info('P{0}     : {1:.4e}'.format(i, v))

#===========================================================================

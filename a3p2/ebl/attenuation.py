#===========================================================================
# Imports

#import math

import numpy as np
import scipy.integrate
import scipy.interpolate
import pyfits

# DEBUG
import matplotlib.pyplot as plt

from a3p2.constants import *

#===========================================================================
# Functions & classes

#---------------------------------------------------------------------------
def convert_ebl_nFnu2n(lambda_micron, nuFnu_nWm2sr1) :
    # Convert EBL data to log10(energy [ J ]) vs log10(number density [ m^-3 ])
    ebl_n_x = 1.239841 / lambda_micron * misc.eV_J
    ebl_n_y = nuFnu_nWm2sr1 * 1E-9 * 4. * np.pi / si.c * ebl_n_x ** -2.
    # Revert axis (lists)
    return (ebl_n_x[::-1], ebl_n_y[::-1])

#---------------------------------------------------------------------------
class EBLModelASCII(object) :

    def __init__(self, filename, ebl_evo_f=0., skiprows=0) :

        self.ebl_evo_f =  ebl_evo_f

        # Read EBL data from ASCII file
        # Format is two column ASCII with
        #  (1) EBL wavelength [ microns ]
        #  (2) EBL intensity [ nW m^-1 sr^-1 ]
        tbdata = np.loadtxt(filename, skiprows=skiprows)
        self.lambda_micron = tbdata[:, 0]
        self.nuFnu_nWm2sr1 = tbdata[:, 1]

        ## Convert EBL data to log10(energy [ J ]) vs log10(number density [ m^-3 ])
        ebl_n_x, ebl_n_y = convert_ebl_nFnu2n(self.lambda_micron, self.nuFnu_nWm2sr1)
        ebl_n_x, ebl_n_y = np.log10(ebl_n_x), np.log10(ebl_n_y)

        # Set EBL model limits
        self.lambda_min_micron, self.lambda_max_micron = min(self.lambda_micron), max(self.lambda_micron)
        self.e_min_J, self.e_max_J = 10.**min(ebl_n_x), 10.**max(ebl_n_x)

        # Fill lookup tables (1d linear interpolation with splines)
        # Sanity
        ebl_n_x[np.invert(np.isfinite(ebl_n_x))] = 0.
        ebl_n_y[np.invert(np.isfinite(ebl_n_y))] = 0.
        self.lookup_n = scipy.interpolate.UnivariateSpline(ebl_n_x, ebl_n_y, k=1, s=0)
        self.lookup_nuFnu = scipy.interpolate.UnivariateSpline(self.lambda_micron, self.nuFnu_nWm2sr1, k=1, s=0)

    def get_n(self, e_J, z=0.) :
        """
        Returns EBL photon number density for a given energy [eV] and redshift
        from a spline interpolation of the data in loglog.
        """
        if not self.has_data :
            return 0.
        if e_J > self.e_max_J or e_J < self.e_min_J :
            return 0.
        else :
            return 10. ** self.lookup_n(np.log10(e_J / (z + 1.))) * (z + 1.) ** (2. - self.ebl_evo_f)# WHY 2. ???

    def get_n_a(self, e_J, z=0.) :
        """
        Returns EBL photon number density for a given energy [eV] and redshift
        from a spline interpolation of the data in loglog.
        """
        return 10. ** self.lookup_n(np.log10(e_J / (z + 1.))) * (z + 1.) ** (2. - self.ebl_evo_f)# WHY 2. ???

#---------------------------------------------------------------------------
class EBLModelBSpline(object) :

    def __init__(self, ebl_lambda_micron, ebl_nuFnu_nWm2sr1, ebl_evo_f=0., useknot=-1, steps_per_dec = 100.) :

        # should use super sampling, i.e. sample spline or base spline and create UnivariateSpline first
        # this speeds up the calculation tremendously
        import a3p2.tools.bspline

        self.ebl_evo_f =  ebl_evo_f
        self.lambda_micron = ebl_lambda_micron
        self.nuFnu_nWm2sr1 = ebl_nuFnu_nWm2sr1

        # Create b-spline for the data
        # B-spline is of order 2 and uses log10(wavelength) vs ebl energy density in nW m^-2 sr^-1
        self.bspline = a3p2.tools.bspline.data_to_bspline1d(np.log10(self.lambda_micron),
                                                            self.nuFnu_nWm2sr1, 2)

        # Convert EBL data to energy [ J ] vs number density [ m^-3 ]
        # This is just used to calculate e_min and e_max - me being lazy
        ebl_n_x, ebl_n_y = convert_ebl_nFnu2n(self.lambda_micron, self.nuFnu_nWm2sr1)

        # Set EBL model limits
        self.lambda_min_micron, self.lambda_max_micron = min(self.lambda_micron), max(self.lambda_micron)
        self.e_min_J, self.e_max_J = min(ebl_n_x), max(ebl_n_x)

        self._useknot = useknot
        self._steps_per_dec = steps_per_dec
        self.sample_bspline()

        #self.get_n_a = np.vectorize(self.get_n)

    def get_useknot(self) :
        return self._useknot
    def set_useknot(self, v) :
        self._useknot = v
        self.sample_bspline()
    useknot = property(get_useknot, set_useknot)
    
    def get_steps_per_dec(self) :
        return self._steps_per_dec
    def set_steps_per_dec(self, steps_per_dec) :
        self._steps_per_dec = steps_per_dec
        self.sample_bspline()
    steps_per_dec = property(get_steps_per_dec, set_steps_per_dec)

    def sample_bspline(self) :
        # Sample spline for EBL number density lookup
        l_tmp = np.linspace(np.log10(self.lambda_min_micron), np.log10(self.lambda_max_micron),
                            np.abs((np.log10(self.lambda_max_micron)
                                    - np.log10(self.lambda_min_micron)) * self._steps_per_dec))
        f = self.bspline
        ebl_x, ebl_y, ebl_n_x, ebl_n_y = None, None, None, None
        if self._useknot == -1 :
            ebl_x, ebl_y = 10. ** l_tmp, np.vectorize(self.bspline.eval)(l_tmp, 2)
        else :
            ebl_x, ebl_y = 10. ** l_tmp, np.vectorize(self.bspline.eval_base_spline)(l_tmp, self._useknot, 2)
            m = ebl_y > 0.
            ebl_x, ebl_y = ebl_x[m], ebl_y[m]
        # Create lookup splines (linear interpolation)
        self.lookup_nuFnu = scipy.interpolate.UnivariateSpline(ebl_x, ebl_y, k=1, s=0)
        ebl_n_x, ebl_n_y = convert_ebl_nFnu2n(ebl_x, ebl_y)
        self.lookup_n = scipy.interpolate.UnivariateSpline(np.log10(ebl_n_x), np.log10(ebl_n_y), k=1, s=0)
    
    def get_n(self, e_J, z=0.) :
        """
        Returns EBL photon number density for a given energy [eV] and redshift
        from a spline interpolation of the data in loglog.
        """
        if self.bspline == None :
            return 0.
        if e_J > self.e_max_J or e_J < self.e_min_J :
            return 0.
        # log10(wavelength [ microns ])
        e_J = e_J / (z + 1.)
        l = np.log10(1.239841 / e_J  / misc.J_eV)
        # [ nW m^-2 sr^-1 ] -> [ m^-3 ]
        c1 = 1E-9 * 4. * np.pi / si.c * e_J ** -2.
        if self.basespline == -1 :
            return self.bspline.eval(l, 2) * c1 * (z + 1.) ** (2. - self.ebl_evo_f)
            # DEBUG DEBUG DEBUG
            #return 10. ** self.bspline.eval(np.log10(e_J / (z + 1.)), 2) * (z + 1.) ** (2. - self.ebl_evo_f)
        else :
            return self.bspline.eval_base_spline(l, self.basespline, 2) * c1 * (z + 1.) ** (2. - self.ebl_evo_f)

    def get_n_a(self, e_J, z=0.) :
        """
        Returns EBL photon number density for a given energy [eV] and redshift
        from a spline interpolation of the data in loglog.
        """
        return 10. ** self.lookup_n(np.log10(e_J / (z + 1.))) * (z + 1.) ** (2. - self.ebl_evo_f)

#---------------------------------------------------------------------------
class EBLModelMR2012(object) :

    def __init__(self, filename) :

        self.ebl_evo_f = 0.
    
        f = pyfits.open(filename)
        dx = f['LOGNU'].data.field('LOGNU') + np.log10(si.h)
        dy = f['EBLZ'].data.field('Z')
        #dz = np.log10(f['EBL'].data * 1E-9 * 4. * np.pi / si.c * (10. ** dx[:,np.newaxis]) ** -2.)
        dz = np.log10(f['EBL'].data) - 9. + np.log10(4. * np.pi / si.c) - 2. * dx[:,np.newaxis]
        
        self.e_min_J, self.e_max_J = 10. ** dx[0], 10. ** dx[-1]

        dz[np.invert(np.isfinite(dz))] = 0.
        self.lookup_n = scipy.interpolate.RectBivariateSpline(x=dx, y=dy, z=dz, kx=1, ky=1)

        self.name = f[0].header['NAME']
        f.close()

    def get_n(self, e_J, z=0.) :
        return get_n_a(e_J, z)

    def get_n_a(self, e_J, z=0.) :
        """
        Returns EBL photon number density for a given energy [eV] and redshift
        from a spline interpolation of the data in loglog.
        """
        return 10. ** self.lookup_n.ev(np.log10(e_J), z) * (z + 1.) ** 3.

#---------------------------------------------------------------------------
def calc_ebl_attenuation(z_end, E_TeV, ebl, h=0.7, W_m=0.3, W_l=0.7) :
    """
    Calculates the optical depth for VHE gamma-rays for a given redshift, energy, and EBL density.
    """

    E_J = E_TeV * 1e12 * misc.eV_J
    LOG10 = np.log(10)

    int_steps_z, int_steps_mu, int_steps_log10e = 21, 21, 501

    z_arr = np.linspace(0., z_end, int_steps_z, endpoint=True)
    mu_arr = np.linspace(-1. + 1E-6, 1. - 1E-6, int_steps_mu, endpoint=True)
    log10_e_arr = np.zeros(int_steps_log10e)
    z_int_arr = np.zeros(int_steps_z)
    mu_int_arr = np.zeros(int_steps_mu)
    log10_e_int_arr = np.zeros(int_steps_log10e)

    # Assign functions/constants to variables to speed up loop
    int_func = scipy.integrate.simps
    #int_func = integrate.trapz
    sqrt, log, log10 = np.sqrt, np.log, np.log10
    linspace = np.linspace
    get_n = ebl.get_n
    get_n_a = ebl.get_n_a
    c, me = si.c, si.me

    for i, z in enumerate(z_arr) :
        for j, mu in enumerate(mu_arr) :
            e_thr = 2. * (c * c * me) ** 2. / E_J / (1. - mu) / (z + 1.)
            log10_e_thr = log10(e_thr)
            log10_e_arr = linspace(log10_e_thr + 1E-8, log10_e_thr + 5., int_steps_log10e,
                                   endpoint=True)

            #for k, log10_e in enumerate(log10_e_arr) :
            #    e = 10. ** log10_e
            #    b = sqrt(1. - e_thr / e)
            #    bb = b * b
            #    r = (1. - bb) * (2. * b * (bb - 2.) + (3. - bb * bb) * log((1. + b) / (1. - b)))
            #    log10_e_int_arr[k] = r * get_n(e, z) * e * LOG10

            # numpy with mask - fastest
            log10_e_int_arr = np.zeros(int_steps_log10e)
            e = 10. ** log10_e_arr
            m = (e < ebl.e_max_J) * (e > ebl.e_min_J)
            if np.any(m) :
                b = np.sqrt(1. - e_thr / e[m])
                bb = b * b
                r = (1. - bb) * (2. * b * (bb - 2.) + (3. - bb * bb) * np.log((1. + b) / (1. - b)))
                #log10_e_int_arr[m] = r * get_n_a(e[m], z) * e[m] * LOG10
                log10_e_int_arr[m] = r * get_n_a(e[m], z * np.ones(np.sum(m))) * e[m] * LOG10

            # numpy
            #b = np.sqrt(1. - e_thr / e)
            #bb = b * b
            #r = (1. - bb) * (2. * b * (bb - 2.) + (3. - bb * bb) * np.log((1. + b) / (1. - b)))
            #log10_e_int_arr = r * get_n_a(e, z) * e * LOG10

            mu_int_arr[j] = int_func(log10_e_int_arr, log10_e_arr) * (1. - mu) / 2.

            # Sanity old
            #if math.isnan(mu_int_arr[j]) :
            #    mu_int_arr[j] = 0.
        # Sanity new
        mu_int_arr[np.invert(np.isfinite(mu_int_arr))] = 0.
        zp1 = z + 1.
        cos = 1E9 * astro.yr_s / (astro.WMAP3_H0100 * h) / (zp1 * sqrt(zp1 * zp1 * (1. + W_m * z) - z * (2. + z) * W_l));
        z_int_arr[i] = int_func(mu_int_arr, mu_arr) * cos
    # Remember: si.c is from dz/dt !
    return 3. / 16. * si.tcs * si.c * int_func(z_int_arr, z_arr)

#---------------------------------------------------------------------------
def calc_ebl_attenuation2(z_end, E_TeV, ebl, h=0.7, W_m=0.3, W_l=0.7,
                          int_steps_z=21, int_steps_mu=21, int_steps_log10e=501) :
    """
    Calculates the optical depth for VHE gamma-rays for a given redshift, energy, and EBL density.
    """

    E_J = E_TeV * 1e12 * misc.eV_J
    LOG10 = np.log(10)

    # Assign functions/constants to variables to speed up loop
    int_func = scipy.integrate.simps
    #int_func = integrate.trapz
    sqrt, log, log10 = np.sqrt, np.log, np.log10
    linspace = np.linspace
    get_n = ebl.get_n
    get_n_a = ebl.get_n_a
    c, me = si.c, si.me

    # Prepare cubes
    cube = np.ones([int_steps_z, int_steps_mu, int_steps_log10e])
    zcube = cube * np.linspace(0., z_end, int_steps_z)[:, np.newaxis, np.newaxis]
    mucube = cube * np.linspace(-1. + 1E-6, 1. - 1E-6, int_steps_mu)[np.newaxis, :, np.newaxis]
    intcube = cube * np.linspace(0., 1., int_steps_log10e)

    lethrcube = np.log10(2. * (c * c * me) ** 2. / E_J / (1. - mucube) / (zcube + 1.))
    lecube = lethrcube + 1E-8 + intcube * 5.

    m = (10. ** lecube < ebl.e_max_J) * (10. ** lecube > ebl.e_min_J)

    # Are we in range of the EBL model?
    if np.sum(m) == 0 :
        return 0.

    intcube[np.invert(m)] = 1E-42

    # Cross section & EBL density
    bcube = np.sqrt(1. - 10. ** lethrcube[m] / 10. ** lecube[m])
    bbcube = bcube * bcube
    r = (1. - bbcube) * (2. * bcube * (bbcube - 2.) + (3. - bbcube * bbcube) * np.log((1. + bcube) / (1. - bcube)))
    intcube[m] = r * get_n_a(10. ** lecube[m].flatten(), zcube[m].flatten()) * 10. ** lecube[m] * LOG10

    # Sanity checks & cleanup
    intcube[np.invert(np.isfinite(intcube)) * np.invert(m)] = 1E-42

    # Integrate over EBL SED
    intcube = int_func(intcube, lecube)

    # Sanity checks & cleanup
    intcube[np.invert(np.isfinite(intcube))] = 1E-42

    # Integrate over cos(theta)
    intcube = int_func(intcube * (1. - mucube[:, :, 0]) / 2., mucube[:, :, 0])

    # Sanity checks & cleanup
    intcube[np.invert(np.isfinite(intcube))] = 1E-42

    # Cosmology
    z = zcube[:, 0, 0]
    zp1 = z + 1.
    cos = 1E9 * astro.yr_s / (astro.WMAP3_H0100 * h) / (zp1 * np.sqrt(zp1 * zp1 * (1. + W_m * z) - z * (2. + z) * W_l));

    # Integrate over redshift
    return 3. / 16. * si.tcs * si.c * int_func(intcube * cos, zcube[:, 0, 0])
    
#===========================================================================

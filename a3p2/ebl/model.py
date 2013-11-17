#===========================================================================
# Imports

import sys
import time
import logging
import gzip

import pyfits
import numpy as np
import scipy.interpolate
import scipy.integrate
import matplotlib.pyplot as plt

import a3p2.astro.cosmo
from a3p2.constants import *
from a3p2.tools.config import ConfigBase

#===========================================================================
# Functions & classes

#---------------------------------------------------------------------------
def read_ssp_sb99_fig(datfile) :
    """
    Read simple stellar population model spectra from the original starburst 99 paper:
    http://www.stsci.edu/science/starburst99/docs/table-index.html
    
    [erg s^-1 A^-1], 1E6 M_solar
    """
    dat = np.loadtxt(datfile)
    dx = np.log10(si.c / dat[1:,0][::-1] / 1E-10) # log(frequency/Hz)
    dy = np.log10(dat[0,1:]) # log(Time/yrs)
    
    dz = (dat[1:, 1:][::-1] - 6. + np.log10(1E10 * si.c) - 2. * dx[:,np.newaxis])# log(em [erg/s/Hz/M_solar])
    # Sanity
    dz[np.invert(np.isfinite(dz))] = -43.
    return (dx, dy, dz, dat)

#---------------------------------------------------------------------------
def read_ssp_sb99(datfile) :
    """
    Read simple stellar population model spectra from starburst 99 output:
    http://www.stsci.edu/science/starburst99/
    
    [erg s^-1 A^-1], 1E6 M_solar
    """
    d = np.loadtxt(datfile, skiprows=6)
    # Get unique time steps
    tt = d[:,0].tolist()
    tt.sort()
    unique = [x for i, x in enumerate(tt) if not i or x != tt[i-1]]
    t = np.asarray(unique)
    # Read spectra, store in data array
    l, dd, first = None, None, True
    for ts in t :
        dt = d[d[:,0] == ts, :]
        if first :
            l = dt[:,1]
            dd = dt[:, 3][:, np.newaxis]
            first = False
        else :
            dd = np.hstack([dd, dt[:, 3][:, np.newaxis]])
    dx = np.log10(si.c / l[::-1] / 1E-10) # log(frequency/Hz)
    dy = np.log10(t) # log(Time/yrs)
    dz = (dd[::-1] - 6. + np.log10(1E10 * si.c) - 2. * dx[:,np.newaxis])# log(em [erg/s/Hz/M_solar])
    # Sanity
    dz[np.invert(np.isfinite(dz))] = -43.
    return (dx, dy, dz, dd)

# BC T. Kneiske et al. 2002
#datfile = a3p2...datapath + '/ssp_salp_z02_bc_all.dat' # [L_solar M_solar^-1 A^-1]
#dat = np.loadtxt(datfile)
#dx = np.log10(si.c / dat[1:,0][::-1] / 1E-10) # log(frequency/Hz)
#dy = np.log10(dat[0,1:]) # log(Time/yrs)
#dz = (np.log10(dat[1:, 1:][::-1] * astro.L_SOLAR * 1E7) + np.log10(1E10 * si.c) - 2. * dx[:,np.newaxis])# log(em [erg/s/Hz/M_solar])

#---------------------------------------------------------------------------
def read_ssp_bc2003(datfile) :
    """
    Read simple stellar population model spectra from Bruzual & Charlot 2003:
    http://www2.iap.fr/users/charlot/bc2003/
    [erg s^-1 A^-1 M_solar^-1]    
    """
    f = None
    if datfile[-3:] == '.gz' :
        f = gzip.open(datfile)
    else :
        f = open(datfile)
    s = f.read()
    f.close()
    s = s.split(' 1221 ')
    l, dat, first, n = None, None, True, 0
    t = np.array([float(x) for x in s[0].split()[1:222]])
    for a in s[1:-1] :
        a = np.array([float(x) for x in a.split()])
        if first :
            first = False
            l = a
        else :
            if n == 0 :
                dat = a[:-53, np.newaxis]
            else :
                dat = np.hstack([dat, a[:-53, np.newaxis]])
            n += 1
    dat = np.hstack([dat, np.array([float(x) for x in s[-1].split(' 221 ')[0].split()[:-53]])[:, np.newaxis]])
    dx = np.log10(si.c / l[::-1] / 1E-10) # log(frequency/Hz)
    if t[0] < 1. :
        t[0] = 1.
    dy = np.log10(t) # log(Time/yrs)
    dz = (np.log10(dat[::-1] * astro.L_SOLAR * 1E7) + np.log10(1E10 * si.c) - 2. * dx[:,np.newaxis])# log(em [erg/s/Hz/M_solar])
    # Sanity
    dz[np.invert(np.isfinite(dz))] = -43.
    return (dx, dy, dz, dat)

#---------------------------------------------------------------------------
def read_ssp_tum(datfile) :
    """
    Read simple stellar population model spectra from Tumlinson 2006.
    [erg s^-1 A^-1 M_solar^-1]    
    """
    dat = np.loadtxt(datfile, skiprows=2)
    dx = np.log10(si.c / dat[1:,0][::-1] / 1E-10) # Wavelength [Angstroem] -> log(frequency/Hz)
    dy = dat[0,1:] # log(Time/yrs)
    dz = (dat[1:, 1:][::-1] + np.log10(1E10 * si.c) - 2. * dx[:,np.newaxis])# log(em/L_solar/Hz)
    # Sanity
    dz[np.invert(np.isfinite(dz))] = -43.
    return (dx, dy, dz, dat)

#---------------------------------------------------------------------------
def calc_emissivity(
    dat_intp,# SSP spectra lookup
    ssp_t, ssp_nu, # log10 SSP time steps/frequency steps
    sfr=lambda zp1: 1.,# star formation rate function f(z+1) !!!!
    lmin=0.08, lmax=500., # Emissivity wavelength min/max [micrometer]
    lfsteps = 201, # Emissivity log10(frequency) steps
    zmin=1E-5, zmax=5., zintmax=5., # Emmissivity redshift minimum/maximum, integration redshift maximum
    zsteps=201,# Emissivity redshift steps
    ltintsteps=201,# log10(time) integration step
    cosmo=[.7, .3, .7], # Cosmology
    **name
    ) :

    z, t = a3p2.astro.cosmo.z_to_t_table(zmin=1E-6, cosmo=cosmo)
    z2t = scipy.interpolate.UnivariateSpline(np.log10(z), np.log10(t) + 9., s=0, k=1) # [log10(yrs)]
    t2z = scipy.interpolate.UnivariateSpline(np.log10(t) + 9., np.log10(z), s=0, k=1) # [log10(yrs)]

    # Setup frequency and redshift ranges
    fmin, fmax = si.c / lmax / 1E-6, si.c / lmin / 1E-6# Frequency min/max [Hz]
    lfmin, lfmax = np.log10([fmin, fmax])
    ltintmax = z2t(np.log10(zintmax))

    logging.info('Lamda range            : {0} - {1} [micrometer]'.format(lmin, lmax))
    logging.info('Frequency range        : {0:.2E} - {1:.2E} [Hz] in {2} log steps'.format(fmin, fmax, lfsteps))
    logging.info('Emissivity z range     : {0} - {1} in {2} steps'.format(zmin, zmax, zsteps))
    logging.info('Emissivity integration : z_max = {0}, ltintsteps = {1}'.format(zintmax, ltintsteps))

    # Setup data cubes
    cube = np.ones([lfsteps, zsteps, ltintsteps])

    lfcube = cube * np.linspace(lfmin, lfmax, lfsteps)[:, np.newaxis, np.newaxis]
    zcube = cube * np.linspace(zmin, zmax, zsteps)[np.newaxis, :, np.newaxis]
    ltintcube = cube * np.linspace(0., 1., ltintsteps)

    logging.info('Calculating t cubes ..')
    ltstartcube = z2t(np.log10(zcube.flatten())).reshape([lfsteps, zsteps, ltintsteps])
    ltdt = np.log10(10. ** ltintmax - 10. ** ltstartcube)
    # limit integration range to spline range
    # Note: for values exceeding the spline boundaries the last value is returned
    ltdt[ltdt > ssp_t[-1]] = ssp_t[-1]
    ltoffset = ssp_t[0]
    logging.debug('ltoffset [log(yr)] = {0}'.format(ltoffset))
    ltintcube = (ltdt - ltoffset) * ltintcube + ltoffset

    # Calculate integration values
    intcube = cube * 1E-43
    # limit integration range to spline range
    s = (lfcube >= ssp_nu[0]) * (lfcube <= ssp_nu[-1])
    logging.debug('Sampling source ..')
    intcube[s] = 10. ** dat_intp.ev(lfcube[s], ltintcube[s]) * 10. ** ltintcube[s] * np.log(10.) * sfr(10. ** t2z(np.log10(10. ** ltstartcube[s] + 10. ** ltintcube[s]).flatten()) + 1.)
    # !!!!! sfr(z+1) !!!!

    # Free memory
    tcube, tstartcube, cube, ltstartcube, ltdt = None, None, None, None, None

    # Emissivity [erg s^-1 Hz^-1 Mpc^-3]
    logging.debug('Integrating ..')
    em = scipy.integrate.simps(intcube, ltintcube) # [erg s^-1 Hz^-1 Mpc^-3]
    lem = np.log10(em)
    lem[np.invert(np.isfinite(lem))] = -43.
    em_intp = scipy.interpolate.RectBivariateSpline(x=np.linspace(lfmin, lfmax, lfsteps),
                                                    y=np.linspace(zmin, zmax, zsteps), z=lem, kx=1, ky=1)
    # Free memory
    lfcube, zcube, intcube, ltintcube, s, lem = None, None, None, None, None, None

    return (em, em_intp)

#---------------------------------------------------------------------------
def calc_ebl(
    em_intp, lmin=0.08, lmax=500., zmax=15., lfsteps=201, zsteps=201,
    eblzmin=0., eblzmax=5., eblzsteps=201, cosmo=[.7, .3, .7], **name
    ) :
    eblzintsteps = zsteps

    fmin, fmax = si.c / lmax / 1E-6, si.c / lmin / 1E-6# Frequency min/max [Hz]
    lfmin, lfmax = np.log10([fmin, fmax])

    cube = np.ones([lfsteps, eblzsteps, eblzintsteps])

    lfcube = cube * np.linspace(lfmin, lfmax, lfsteps)[:, np.newaxis, np.newaxis]
    eblzcube = cube * np.linspace(eblzmin, eblzmax, eblzsteps)[np.newaxis, :, np.newaxis]
    eblzintcube = cube * np.linspace(0., 1., eblzintsteps)

    eblzintcube = eblzcube + eblzintcube * (zmax - eblzcube)

    # Calculate integration values
    eblintcube = 10. ** em_intp.ev((lfcube + np.log10((1. + eblzintcube) / (1. + eblzcube))).flatten(),
                                   eblzintcube.flatten()).reshape([lfsteps, eblzsteps, eblzintsteps])
    eblintcube *= a3p2.astro.cosmo.cosmo_term1(*cosmo)(eblzintcube) * 1E9 # Gyr -> yrs
    #yr -> s, Mpc^-3 -> m^-3, erg -> nJ,
    eblintcube *= 1.068885e-58
    eblintcube *= 10. ** lfcube * si.c / 4. / np.pi

    # EBL
    ebl = scipy.integrate.simps(eblintcube, eblzintcube)
    lebl = np.log10(ebl)
    lebl[np.isnan(lebl)] = -43.
    lebl[np.invert(np.isfinite(lebl))] = -43.
    ebl_intp = scipy.interpolate.RectBivariateSpline(x=np.linspace(lfmin, lfmax, lfsteps),
                                                     y=np.linspace(eblzmin, eblzmax, eblzsteps),
                                                     z=lebl, kx=1, ky=1)
    # Free memory
    lfcube, eblzcube, eblzintcube, eblintcube, cube, lebl = None, None, None, None, None, None

    return (ebl, ebl_intp)

#---------------------------------------------------------------------------
def nebula_emission(l, calcff=True, calcfb=True, calc2ph=True, calclya=True) :
    """
    Stellar nebula emission following Fernandez & Komatsu 2006.
    """

    nu = si.c / (l * 1E-10) # [Hz] Frequency

    em = np.zeros(len(l))

    nuLyAlpha = si.c / 1215.67E-10 # [Angstroem]

    # Free-bound continuum emission
    # Gaunt factors
    Tg = 20000. # [K] Temperature
    fbconsta = si.me * si.c ** 2. * si.alpha ** 2. / 2. / si.k / Tg
    nbase = np.arange(2., 10.)[:, np.newaxis] + np.arange(0., 50.)[np.newaxis, :]
    fbfact = fbconsta * 1.05 * np.sum(np.exp(fbconsta / nbase ** 2.) / nbase ** 3., axis=1)

    # Free-free continuum emission
    # Fernandez & Komatsu 2006, Eq. 12
    if calcff :
        em += 1.1

    # Free-bound continuum emission
    # Fernandez & Komatsu 2006, Eq. 12
    if calcfb :
        llim = [3646., 8201., 14580., 22782., 32810.] # [Angstroem]
        limold = 0.
        for lim, fbf in zip(llim, fbfact) :
            em = np.where((l > limold) * (l <= lim), em + fbf, em)
            limold = lim
        em = np.where(l > llim[-1], em + fbfact[5], em)

    em *= 3.32E22 / 1e49 * np.exp(-si.h * nu / si.k / Tg);
    #      em *= 3.32E22*n0/1e49*std::exp(-SI_h*nu/SI_k/Tg);
    #      // erg s^-1 Hz^-1

    # 2 photon emission
    if calc2ph :
        y = nu / nuLyAlpha
        em2p = 2. * 6.626 * y * (1.307 - 2.627 * (y - .5) ** 2. + 2.563 * (y - .5) ** 4. - 51.69 * (y - .5) ** 6.)
        em2p = np.where(em2p < 0., 0., em2p)
        em += em2p * 1E-27 * .33
        # erg s^-1 Hz^-1

    # Ly alpha emission
    if calclya :
        idx=(np.abs(nu - nuLyAlpha)).argmin()
        if idx != 0 and idx != len(nu) - 1 :
            em[idx] += 0.64 * si.h * misc.J_erg
            # erg s^-1 Hz^-1

    return em

#---------------------------------------------------------------------------
def apply_dust(d, ebv=.15, ebv2=None, ir_fac=3e9, ir_wave_start=5.5) :
    dx, dy, dz, dat = d
    # Resample data to encompass larger wavelength range
    logging.info('Extending wavelength range of SSP models for IR emission')
    dat_intp = scipy.interpolate.RectBivariateSpline(x=dx, y=dy, z=dz, kx=1, ky=1)
    dx_old = dx
    ir_lnu_min, ir_lnu_max = 11., 17.
    ir_lnu_nsteps = int((ir_lnu_max - ir_lnu_min) / ((dx_old[-1] - dx_old[0]) / len(dx_old)))
    logging.info('ir_lnu_min={0}, ir_lnu_max={1}, ir_lnu_nsteps={2}'.format(ir_lnu_min, ir_lnu_max, ir_lnu_nsteps))
    dx = np.linspace(ir_lnu_min, ir_lnu_max, ir_lnu_nsteps)
    xx = np.ones([len(dx), len(dy)]) * dx[:, np.newaxis]
    yy = np.ones([len(dx), len(dy)]) * dy[np.newaxis, :]
    dz = dat_intp.ev(xx.flatten(), yy.flatten()).reshape([len(dx), len(dy)])
    m = (dx < dx_old[0]) + (dx > dx_old[-1])
    dz[m] = -43.
    dat_intp = scipy.interpolate.RectBivariateSpline(x=dx, y=dy, z=dz, kx=1, ky=1)

    def att_kn2002(l) :
        return -.4 *  .68 * 3.2 * ((1. / l) - .35)

    def att_fm2007(l, x0=4.592, g=.922, c1=-.175, c2=.807, c3=2.991, c4=.319, c5=6.097,
                   O1=2.055, O2=1.322, O3=0., RV=3.001, kir=1.057) :
        x = 1. / l
        dr = lambda x_, x0_, g_ : x_ * 2. / ((x_ ** 2.  - x0_ ** 2.) ** 2. + x_ ** 2. * g_ ** 2.)
        def fk(x) :
            k = c1 + c2 * x + c3 * dr(x, x0, g)
            m = x > c5
            if np.any(m) :
                k[m] = c1 + c2 * x + c3 * dr(x, x0, g) + c4 * (x - c5) ** 2.
            return k
        k = fk(x)
        m = x < 1. / .27
        if np.any(m) :
            lU1, lU2 = .26, .27
            lO1, lO2, lO3 = .33, .4, .553
            xIR = [0.0001, .25, .5, .75, 1.]
            xI = np.array([0.0001, .25, .5, .75, 1., 1. / .533, 1. / .4, 1. / .33, 1. / .27, 1. / .26])
            yI = np.zeros(len(xI))
            for vi, v in enumerate(xI[:5]) :
                yI[vi] = kir * ((1. / v) ** -1.84) - RV
            yI[5], yI[6], yI[7] = O3, O2, O1
            yI[8], yI[9] = fk(xI[-2:])
            intp = scipy.interpolate.UnivariateSpline(x=xI, y=yI, s=0, k=2)
            k[m] = intp(x[m])
        return -.4 * (k + RV)

        
    # Calculate dust attenuation
    dz_att = np.copy(dz)
    dx_l = (si.c / (10. ** dx) * 1E6)
    if ebv2 :
        logging.info('Using two step dust emission a la KN02: E(B-V)_1 = {0}, E(B-V)_2 = {1}'.format(ebv, ebv2))
        m = dy <= np.log10(3E8)
        dz_att[:, m] += (ebv * att_fm2007(dx_l))[:, np.newaxis]
        #dz_att[:, m] += (ebv * att_kn2002(dx_l))[:, np.newaxis]
        #dz_att[:, m] += (-.4 *  0.68 * ebv * 3.2 * (1. / (si.c / (10. ** dx) * 1E6) - .35))[:, np.newaxis]
        m = dy > np.log10(3E8)
        dz_att[:, m] += (ebv * att_fm2007(dx_l))[:, np.newaxis]
        #dz_att[:, m] += (ebv2 * att_kn2002(dx_l))[:, np.newaxis]
        #dz_att[:, m] += (-.4 *  0.68 * ebv2 * 3.2 * (1. / (si.c / (10. ** dx) * 1E6) - .35))[:, np.newaxis]
    else :
        dz_att += (ebv * att_fm2007(dx_l))[:, np.newaxis]
        #dz_att += (ebv * att_kn2002(dx_l))[:, np.newaxis]
        #dz_att += (-.4 *  0.68 * ebv * 3.2 * (1. / (si.c / (10. ** dx) * 1E6) - .35))[:, np.newaxis]


    m = (dx < np.log10(13.6 * misc.eV_J / si.h)) * (dx > np.log10(si.c / ir_wave_start / 1E-6))
    int1 = scipy.integrate.simps(
        np.transpose(10. ** dz_att[m, :] * (10. ** dx[m])[:, np.newaxis] * np.log(10.)),
        np.transpose(np.ones(dz_att[m, :].shape) * dx[m][:, np.newaxis])
        )
    # DEBUG
    # Double UV
    #m =  (dx > np.log10(si.c / ir_wave_start / 1E-6))
    # END DEBUG
    int2 = scipy.integrate.simps(
        np.transpose(10. ** dz[m, :] * (10. ** dx[m])[:, np.newaxis] * np.log(10.)),
        np.transpose(np.ones(dz[m, :].shape) * dx[m][:, np.newaxis])
        )
    int1 *= misc.erg_J / astro.L_SOLAR
    int2 *= misc.erg_J / astro.L_SOLAR
    l_dust_att = int2 - int1
    #DEBUG
    #l_dust_att[dy < 6.75] *= 2.
    #END DEBUG
    m = dx > np.log10(13.6 * misc.eV_J / si.h)
    l_ion = scipy.integrate.simps(
        np.transpose(10. ** dz[m] * (10. ** dx[m])[:, np.newaxis] * np.log(10.)),
        np.transpose(np.ones(dz[m].shape) * dx[m][:, np.newaxis])
        ) *  misc.erg_J / astro.L_SOLAR
    f_lya = .68
    ir_fac = 1. / ir_fac

    # Read IR SEDs from Chary & Elbaz 2001
    f = pyfits.open(a3p2.datapath + '/ssp/chary_elbaz_2001/chary_elbaz.fits.gz')
    d = f[1].data
    ir_l = d.field('LAMBDA')[0]
    ir_lnu = np.log10(si.c / ir_l[::-1] / 1E-6)
    ir_lf = np.log10(d.field('NULNUINLSUN')[0])
    x, y = ir_lf.shape
    ir_lf = ir_lf[::-1]  - ir_lnu[:, np.newaxis] + np.log10(ir_fac) + np.log10(astro.L_SOLAR * 1E7)              
    ir_lint = np.log10(d.field('LIR')[0]) + np.log10(ir_fac)
    # Calculate integrated IR emission from ir_wave_start to 1E3 micrometer
    m = (ir_lnu < np.log10(si.c / ir_wave_start / 1E-6)) * (ir_lnu > np.log10(si.c / 1E3 / 1E-6))
    ir_lint_alt = np.log10(
        scipy.integrate.simps(np.transpose(10. ** (ir_lf[m] - 12.)* (10. ** ir_lnu[m])[:, np.newaxis] * np.log(10.)),
                              np.transpose(np.ones(ir_lf[m].shape) * ir_lnu[m][:, np.newaxis]))
        * misc.erg_J / astro.L_SOLAR
        ) + 12.

    ir_lint = ir_lint_alt
    ir_intp = scipy.interpolate.RectBivariateSpline(x=ir_lnu, y=ir_lint, z=ir_lf, kx=1, ky=1)

    ir_lnu_cube = np.ones(dz.shape) * dx[:, np.newaxis]
    ir_lint_cube = np.ones(dz.shape) *  np.log10((l_dust_att + l_ion * f_lya)[np.newaxis, :])

    ir_dz = ir_intp.ev(ir_lnu_cube.flatten(), ir_lint_cube.flatten()).reshape(dz.shape)

    # Cut off IR emission at wavelength ir_wave_start
    ir_dz[dx > np.log10(si.c / ir_wave_start / 1E-6)] = -43.

    # rescale IR emission for out of range spectra
    m = ir_lint_cube > ir_lint[-1]
    if np.any(m) :
        logging.info('IR emission: scaling up dust emission for {0} spectra'.format(np.sum(m) / len(dx)))
        ir_dz[m] += (ir_lint_cube[m] - ir_lint[-1])
    m = ir_lint_cube < ir_lint[0]
    if np.any(m) :
        logging.info('IR emission: scaling down dust emission for {0} spectra'.format(np.sum(m) / len(dx)))
        ir_dz[m] += (ir_lint_cube[m] - ir_lint[0])

    # !!!!!!!!!!!!!!!!!!!
    # Set emission to infrared + attenuated UV-O-NIR emission
    dz = np.log10(10. ** ir_dz + 10. ** dz_att)

    # DEBUG
    #m = (dx < np.log10(si.c / ir_wave_start / 1E-6)) * (dx > np.log10(si.c / 1E3 / 1E-6))
    #int_debug = scipy.integrate.simps(
    #    np.transpose(10. ** dz[m, :] * (10. ** dx[m])[:, np.newaxis] * np.log(10.)),
    #    np.transpose(np.ones(dz[m, :].shape) * dx[m][:, np.newaxis])
    #    ) * misc.erg_J / astro.L_SOLAR
    #
    #plt.figure()
    #plt.hist((int_debug / (l_dust_att  + l_ion * f_lya)).flatten(), bins=20)
    #plt.figure()
    #plt.plot(dy, (int_debug / (l_dust_att  + l_ion * f_lya)).flatten())
    #plt.figure()
    #plt.semilogx(int_debug, (int_debug / (l_dust_att  + l_ion * f_lya)).flatten())
    #plt.figure()
    #plt.semilogy(dy, l_dust_att / misc.erg_J * astro.L_SOLAR, 'o-')
    #plt.semilogy(dy, int1 / misc.erg_J * astro.L_SOLAR, 'o-')
    #plt.semilogy(dy, int2 / misc.erg_J * astro.L_SOLAR, 'o-')
    #at = np.loadtxt('tmp.txt')
    #plt.semilogy(np.log10(at[:,3]), at[:,0], 'x', color='red')
    #plt.semilogy(np.log10(at[:,3]), at[:,1], 'x', color='green')
    #plt.semilogy(np.log10(at[:,3]), at[:,2], 'x', color='blue')
    #plt.xlabel('log(t)')
    #plt.ylabel('Emissivity (erg/s)')
    #plt.show()
    # DEBUG END

    # Release memory
    dat_intp = None
    dz_att = None
    ir_l, ir_lnu, ir_lf, d = None, None, None, None
    ir_lint, ir_intp = None, None
    ir_lnu_cube, ir_lint_cube, ir_dz = None, None, None
    f.close()
    
    return (dx, dy, dz, dat)

#---------------------------------------------------------------------------
def apply_nebula(d, fesc=0.) :

    if fesc < 0. :
        logging.error('fesc < 0. ({0}) : using fesc=0. for nebula emission'.format(fesc))
        fesc = 0.
    elif fesc > 1. :
        logging.warning('fesc > 1. ({0}) : using fesc=1. for nebula emission (i.e. no nebula emission)'.format(fesc))
        return d
    
    dx, dy, dz, dat = d    
    m = dx > np.log10(13.6 * misc.eV_J / si.h)
    n_ion = scipy.integrate.simps(
        np.transpose(10. ** dz[m] *  misc.erg_J * (10. ** dx[m])[:, np.newaxis] * np.log(10.)
                     / ((10. ** dx[m])[:, np.newaxis] * si.h)),
        np.transpose(np.ones(dz[m].shape) * dx[m][:, np.newaxis])
        ) * (1. - fesc)

    dx_neb = np.copy(dx)
    if dx[0] > 10.5 :
        dx_neb = np.hstack((np.arange(10.5, dx[0], 0.1), dx))
        dz_ext = np.ones((len(dx_neb), len(dy)))
        dz_ext[-len(dx):, :] = dz
        dz = dz_ext
    nebem = nebula_emission(si.c / 10. ** dx_neb * 1E10)
    dz_neb = np.log10(n_ion[np.newaxis, :] * nebem[:, np.newaxis])

    dz = np.log10(10. ** dz + 10. ** dz_neb)
    dx = dx_neb
    # !! ssp_nu is used to dermine frequency range for emissivity calculation
    #p['ssp_nu'] = dx

    # Free memory
    dz_neb, n_ion, nebem = None, None, None
    
    return (dx, dy, dz, dat)

#---------------------------------------------------------------------------
CONFIG_EBL_MODEL_DEFAULT = {
    'name' : 'default',
    'sfrf' : 'lambda p, x : np.where(x < p[0], p[1] * (x / p[0]) ** -p[2],  p[1] * (x / p[0]) ** -p[3])',
    'sfrp0' : [2.1, .15, -3.4, 0.],
    'spsfile' : '/Users/mraue/Stuff/work/ebl/eblmodel-ir-new/bruzual_charlot_2003/bc03/models/Padova1994/chabrier/bc2003_lr_m62_chab_ssp.ised_ASCII',
    'spstype' : 'BC2003',
    'lmin' : 0.08,
    'lmax' : 3E3,
    'lfsteps' : 201,
    'zmin' : 1E-5,
    'zmax' : 5.,
    'zintmax' : 5.,
    'zsteps' : 201,
    'ltintsteps' : 201,
    'eblzmin' : 0.,
    'eblzmax' : 5.,
    'eblzsteps' : 101,
    #'ssp_t' : d[1],
    #'ssp_nu' : d[0],
    'cosmo' : [.7, .3, .7],
    #'Ebv' : .14,
    'ebv' : .15,
    'ebv2' : None,
    #'Ebv2' : .03,
    'ir_wave_start' : 5.5,#5.5
    'ir_fac' : 3E9,
    'dust' : True,
    'nebula' : False,
    'fesc' : 0.
    }

#---------------------------------------------------------------------------
class ConfigEBLModel(ConfigBase) :
    def __init__(self, dict_=None, verbose=True) :
        ConfigBase.__init__(self, CONFIG_EBL_MODEL_DEFAULT)
        if dict_ : ConfigBase.init_from_dict(self, dict_, verbose)

#---------------------------------------------------------------------------
class EBLModel(object) :

    def __init__(self, config) :
        self.config = config
        self.data = {}
        spsreader = read_ssp_bc2003
        if config.spstype == 'SB99' :
            spsreader = read_ssp_sb99
        elif config.spstype == 'TUM2006' :
            spsreader = read_ssp_tum
        self.data['sps'] = spsreader(config.spsfile)
        self.set_sfr(config.sfrf, config.sfrp0)
        #self.data['sfr'] = lambda x: eval(config.sfrf)(config.sfrp0, x) # [M_solar Mpc^-3] f(z+1)

    def set_sfr(self, sfrf, sfrfp0) :
        self.data['sfr'] = lambda x: eval(sfrf)(sfrfp0, x) # [M_solar Mpc^-3] f(z+1)

    def calculate(self, resample=True, fesc=0.) :
        c = self.config
        dx, dy, dz, dat = self.data['sps']
        dx, dy, dz, dat = np.copy(dx), np.copy(dy), np.copy(dz), np.copy(dat)
        #DEBUG
        #ssp_select = [1 , 21 , 28 , 42 , 49 , 56 , 63 , 70 , 77 , 84 , 91 , 99 , 105 , 113 , 118 , 123 , 125 , 127 , 132 , 139 , 145 , 150 , 155 , 165 , 166 , 173 , 185 , 193 , 209]
        #ssp_select = [x - 1 for x in ssp_select]
        #dy = dy[ssp_select]
        #dz = dz[:, ssp_select]
        # DEBUG
        #dz[:, 3] *= 1.002
        # END DEBUG

        if c.nebula :
            # Apply nebula emission
            logging.info('Applying nebula emission')
            dx, dy, dz, dat = apply_nebula((dx, dy, dz, dat), fesc)

        if c.dust :
            # Apply dust attenuation and emission
            logging.info('Applying dust attenuation & emission')
            dx, dy, dz, dat = apply_dust((dx, dy, dz, dat), ebv=c.ebv, ebv2=c.ebv2,
                                         ir_fac=c.ir_fac, ir_wave_start=c.ir_wave_start)

        # Remove ionizing emission
        m = dx > np.log10(13.6 * misc.eV_J / si.h)
        if fesc == 0. :
            dz[m, :] = -43.
        else :
            dz[m, :] += np.log10(fesc)

        # SPS
        dat_intp = scipy.interpolate.RectBivariateSpline(x=dx, y=dy, z=dz, kx=1, ky=1)

        if resample :
            # Resample data in coarser steps in log(Freq.)
            logging.info('Resampling data in frequency')
            dx = np.linspace(dx[0], dx[-1], 600)
            xx = np.ones([len(dx), len(dy)]) * dx[:, np.newaxis]
            yy = np.ones([len(dx), len(dy)]) * dy[np.newaxis, :]
            dz = dat_intp.ev(xx.flatten(), yy.flatten()).reshape([len(dx), len(dy)])
            dat_intp = scipy.interpolate.RectBivariateSpline(x=dx, y=dy, z=dz, kx=1, ky=1)

        logging.debug('INPUT SSP: lfmin = {0}, lfmax={1}, lfnsteps={2}'.format(dx[0], dx[-1], len(dx)))
        logging.debug('INPUT SSP: ltmin = {0}, ltmax={1}, ltnsteps={2}'.format(dy[0], dy[-1], len(dy)))

        #---------------------------------------------------
        # Calculate emissivity
        logging.info('Calculating emissivity ..')
        cdir = c.__dict__
        cdir['dat_intp'], cdir['ssp_t'], cdir['ssp_nu']= dat_intp, dy, dx
        cdir['sfr'] = self.data['sfr']
        self.data['em'], self.data['em_intp'] = calc_emissivity(**cdir)

        #---------------------------------------------------------------------------
        # Calculate EBL
        logging.info('Calculating EBL ..')
        cdir['em_intp'] = self.data['em_intp']
        self.data['ebl'], self.data['ebl_intp'] = calc_ebl(**cdir)

        # Some cleanup
        dx, dy, dz, dat = None, None, None, None
        
        logging.info('done.')

    def write_fits(self, filename) :
        if len(self.data) == 0 :
            logging.error( 'Could not create HDUs: data has not yet been calculated.')
            return
        lf = np.linspace(
            np.log10(si.c / self.config.lmax / 1E-6),
            np.log10(si.c / self.config.lmin / 1E-6),
            self.config.lfsteps
            )
        hdulist =  pyfits.HDUList([
            pyfits.PrimaryHDU(),
            pyfits.ImageHDU(data=self.data['em'], name='EM'),
            pyfits.ImageHDU(data=self.data['ebl'],name='EBL'),
            pyfits.new_table(pyfits.ColDefs([
                pyfits.Column(name='LOGNU',format='E', array=lf, unit='log10(Hz)')
                ])),
            pyfits.new_table(pyfits.ColDefs([
                pyfits.Column(name='WAVE',format='E', array=si.c / 10. ** lf * 1E6, unit='micrometer')
                ])),
            pyfits.new_table(pyfits.ColDefs([
                pyfits.Column(name='Z',format='E', array=np.linspace(self.config.zmin,
                                                                     self.config.zmax,
                                                                     self.config.zsteps))
                ])),
            pyfits.new_table(pyfits.ColDefs([
                pyfits.Column(name='Z',format='E', array=np.linspace(self.config.eblzmin,
                                                                     self.config.eblzmax,
                                                                     self.config.eblzsteps))
                ]))
            ])
        hdulist[3].header.update('EXTNAME', 'LOGNU')
        hdulist[4].header.update('EXTNAME', 'WAVEMICR')
        hdulist[-2].header.update('EXTNAME', 'EMZ')
        hdulist[-1].header.update('EXTNAME', 'EBLZ')
        self.config.write_to_fits_header(hdulist[0].header,
                                         exclude=['sfr', 'dat_intp', 'ssp_t', 'ssp_nu', 'sfr', 'em_intp'],
                                         rename={'ir_wave_start' : 'irwavest'})
        hdulist.writeto(filename)

#===========================================================================
#===========================================================================

#===========================================================================
# Imports

import numpy as np
import matplotlib.pyplot as plt

#===========================================================================
# Functions & classes

#---------------------------------------------------------------------------
def get_ebl_measurement_collection() :
    ebl_limit_collection = [
        {
            'full_name': 'Bernstein et al. 2002, 2005',
            'id': 'bernstein2002',
            'reference': '',
            'year': 2002,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [0.3000, 0.5550, 0.8140],
                'ebl': [12.0000, 15.0000, 17.9000],
                'ebl_err_low': [7.5000, 7.7000, 8.1000],
                'ebl_err_high': [7.5000, 7.7000, 8.1000]
                }
            }
        ,
        {
            'full_name': 'Brown et al. 2000 (HST/STIS)',
            'id': 'brown2000',
            'reference': 'AJ, 120, 1153-59',
            'year': 2000,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [0.1595],
                'ebl': [14.0000],
                'ebl_err_low': [0.0000],
                'ebl_err_high': [0.0000]
                }
            }
        ,
        {
            'full_name': 'Cambresy et al. 2001 (DIRBE/2MASS)',
            'id': 'cambresy2001',
            'reference': 'ApJ, 555, 2, pp. 563-571.',
            'year': 2001,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [1.2500, 2.2000],
                'ebl': [54.0000, 28.0000],
                'ebl_err_low': [17.0000, 7.0000],
                'ebl_err_high': [17.0000, 7.0000]
                }
            }
        ,
        {
            'full_name': 'Dole et al. 2004 (SPITZER)',
            'id': 'dole2004',
            'reference': 'ApJS, 154',
            'year': 2004,
            'is_lower_limit': 1,
            'is_upper_limit': 0,
            'data': {
                'lambda': [70.0000, 160.0000],
                'ebl': [0.9500, 1.4000],
                'ebl_err_low': [0.0000, 0.0000],
                'ebl_err_high': [0.0000, 0.0000]
                }
            }
        ,
        {
            'full_name': 'Dole et al. 2006 (SPITZER)',
            'id': 'dole2006',
            'reference': 'A&A, 451, 417',
            'year': 2006,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [70.0000, 160.0000],
                'ebl': [7.1000, 13.4000],
                'ebl_err_low': [1.0000, 1.7000],
                'ebl_err_high': [1.0000, 1.7000]
                }
            }
        ,
        {
            'full_name': 'Dube 1979/Leinert 1998',
            'id': 'dube1979',
            'reference': 'Hauser & Dwek 2001',
            'year': 1979,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [0.5115],
                'ebl': [48.0000],
                'ebl_err_low': [0.0000],
                'ebl_err_high': [0.0000]
                }
            }
        ,
        {
            'full_name': 'Dwek & Arendt 1998 (DIRBE)',
            'id': 'dwekArendt1998',
            'reference': '',
            'year': 1998,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [3.5000],
                'ebl': [14.0000],
                'ebl_err_low': [3.0000],
                'ebl_err_high': [3.0000]
                }
            }
        ,
        {
            'full_name': 'Dwek & Arendt UL 1998 (DIRBE)',
            'id': 'dwekArendtUL1998',
            'reference': 'pJ, 508, L9',
            'year': 1998,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [1.2500, 4.9000],
                'ebl': [108.0000, 38.0000],
                'ebl_err_low': [0.0000, 0.0000],
                'ebl_err_high': [0.0000, 0.0000]
                }
            }
        ,
        {
            'full_name': 'Edelstein et al. 2000 (Shuttle/UVX)',
            'id': 'edelstein2000',
            'reference': '',
            'year': 2000,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [0.1000],
                'ebl': [11.0000],
                'ebl_err_low': [0.0000],
                'ebl_err_high': [0.0000]
                }
            }
        ,
        {
            'full_name': 'Elbaz et al. 2002 (ISO)',
            'id': 'elbaz2002',
            'reference': 'A&A, 381',
            'year': 2002,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [15.0000],
                'ebl': [2.4000],
                'ebl_err_low': [0.4000],
                'ebl_err_high': [0.4000]
                }
            }
        ,
        {
            'full_name': 'Fazio et al. 2004 (SPITZER)',
            'id': 'fazio2005',
            'reference': 'Kashlinksy 2005, Physics Report 409 (2005) 361-438',
            'year': 2004,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [3.6000, 4.5000, 5.8000, 8.0000],
                'ebl': [5.2700, 3.9500, 2.7300, 2.4600],
                'ebl_err_low': [1.0200, 0.7700, 0.2200, 0.2100],
                'ebl_err_high': [1.0200, 0.7700, 0.2200, 0.2100]
                }
            }
        ,
        #{
        #    'full_name': 'Fazio et al. 2004 (SPITZER)',
        #    'id': 'fazio2004',
        #    'reference': '',
        #    'year': 2004,
        #    'is_lower_limit': 0,
        #    'is_upper_limit': 0,
        #    'data': {
        #        'lambda': [3.6000, 4.5000, 5.8000, 8.0000],
        #        'ebl': [5.4, 3.5, 3.6, 2.6],
        #        'ebl_err_low': [x * .2 for x in [5.4, 3.5, 3.6, 2.6]],#[.43, .28, .28, .2],
        #        'ebl_err_high':[x * .2 for x in [5.4, 3.5, 3.6, 2.6]]#[.43, .28, .28, .2]
        #        }
        #    }
        #,
        {
            'full_name': 'Finkbeiner et al. 2000 (DIRBE)',
            'id': 'finkbeiner2000',
            'reference': 'ApJ, 544, 81',
            'year': 2000,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [60.0000, 100.0000],
                'ebl': [28.0000, 25.0000],
                'ebl_err_low': [7.0000, 8.0000],
                'ebl_err_high': [7.0000, 8.0000]
                }
            }
        ,
        {
            'full_name': 'Frayer et al. 2006 (SPITZER)',
            'id': 'frayer2006',
            'reference': 'ApJ, 131, 250',
            'year': 2006,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [71.4000],
                'ebl': [7.4000],
                'ebl_err_low': [1.9000],
                'ebl_err_high': [1.9000]
                }
            }
        ,
        {
            'full_name': 'Gorjian et al. 2000 (DIRBE)',
            'id': 'gorjian2000',
            'reference': 'ApJ, 536, 550',
            'year': 2000,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [2.2000, 3.5000],
                'ebl': [22.0000, 11.0000],
                'ebl_err_low': [6.0000, 3.0000],
                'ebl_err_high': [6.0000, 3.0000]
                }
            }
        ,
        {
            'full_name': 'Hauser et al. 1998 (DIRBE/FIRAS)',
            'id': 'hauser1998',
            'reference': 'ApJ, 508, 25',
            'year': 1998,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [140.0000, 140.0000, 240.0000, 240.0000],
                'ebl': [25.0000, 15.0000, 14.0000, 13.0000],
                'ebl_err_low': [7.0000, 6.0000, 3.0000, 2.0000],
                'ebl_err_high': [7.0000, 6.0000, 3.0000, 2.0000]
                }
            }
        ,
        {
            'full_name': 'Hauser et al. UL 1998 (DIRBE/FIRAS)',
            'id': 'hauser1998UL',
            'reference': 'ApJ, 508, 25',
            'year': 1998,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [1.2500, 2.2000, 3.5000, 4.9000, 12.0000, 25.0000, 60.0000, 100.0000],
                'ebl': [75.0000, 39.0000, 23.0000, 41.0000, 470.0000, 500.0000, 75.0000, 34.0000],
                'ebl_err_low': [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
                'ebl_err_high': [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000]
                }
            }
        ,
        {
            'full_name': 'Kashlinsky et al. 1996',
            'id': 'kashlinksy1996',
            'reference': 'ApJ, 470, 681',
            'year': 1996,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [1.2500, 2.2000, 3.5000],
                'ebl': [200.0000, 78.0000, 26.0000],
                'ebl_err_low': [0.0000, 0.0000, 0.0000],
                'ebl_err_high': [0.0000, 0.0000, 0.0000]
                }
            }
        ,
        #{
        #    'full_name': 'Kashlinsky & Odenwald 2000',
        #    'id': 'kashlinskyOdenwald2000',
        #    'reference': 'ApJ, 528, 74',
        #    'year': 2000,
        #    'is_lower_limit': 0,
        #    'is_upper_limit': 1,
        #    'data': {
        #        'lambda': [12.0000, 25.0000, 60.0000, 100.0000],
        #        'ebl': [15.0000, 8.0000, 12.0000, 17.0000],
        #        'ebl_err_low': [0.0000, 0.0000, 0.0000, 0.0000],
        #        'ebl_err_high': [0.0000, 0.0000, 0.0000, 0.0000]
        #        }
        #    }
        #,
        {
            'full_name': 'Lagache et al. 2000 (DIRBE)',
            'id': 'lagache2000',
            'reference': 'AA, 354, 247',
            'year': 2000,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [100.0000],
                'ebl': [23.0000],
                'ebl_err_low': [6.0000],
                'ebl_err_high': [6.0000]
                }
            }
        ,
        {
            'full_name': 'Lagache et al. UL 2000 (DIRBE)',
            'id': 'lagacheUL2000',
            'reference': 'AA, 354, 247',
            'year': 2000,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [140.0000, 240.0000],
                'ebl': [47.0000, 25.0000],
                'ebl_err_low': [0.0000, 0.0000],
                'ebl_err_high': [0.0000, 0.0000]
                }
            }
        ,
        {
            'full_name': 'Levenson & Wright 2008 (SPITZER)',
            'id': 'levenson-wright2008',
            'reference': 'arXiv:0802.1239',
            'year': 2008,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [3.6000, 3.6000],
                'ebl': [9.0000, 7.6000],
                'ebl_err_low': [0.9000, 0.6000],
                'ebl_err_high': [1.7000, 1.0000]
                }
            }
        ,
        {
            'full_name': 'Levenson et al. 2007 (DIRBE/2MASS)',
            'id': 'levenson2007',
            'reference': 'ApJ, 666, 34-44',
            'year': 2007,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [1.2500, 2.2000, 3.5000],
                'ebl': [21.3600, 20.0500, 13.3700],
                'ebl_err_low': [15.1200, 6.1400, 2.8300],
                'ebl_err_high': [15.1200, 6.1400, 2.8300]
                }
            }
        ,
        {
            'full_name': 'Madau & Pozzetti 2000 (HST)',
            'id': 'madauPozzetti2000',
            'reference': 'MNRAS, 312, 2, pp. L9-L15.',
            'year': 2000,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [0.3600, 0.4500, 0.6700, 0.8100, 1.1000, 1.6000, 2.2000],
                'ebl': [2.8700, 4.5700, 6.7400, 8.0400, 9.7000, 9.0000, 7.9000],
                'ebl_err_low': [0.4200, 0.4700, 0.9400, 0.9200, 1.9000, 1.7000, 1.2000],
                'ebl_err_high': [0.5800, 0.7300, 1.2500, 1.6200, 3.0000, 2.6000, 2.0000]
                }
            }
        ,
        {
            'full_name': 'Martin et al. 1991 (Shuttle/UVX)',
            'id': 'martin1991',
            'reference': 'Hauser & Dwek 2001',
            'year': 1991,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [0.1650],
                'ebl': [7.0000],
                'ebl_err_low': [0.0000],
                'ebl_err_high': [0.0000]
                }
            }
        ,
        {
            'full_name': 'Matsumoto et al. 2005 (IRTS)',
            'id': 'matsumoto2005',
            'reference': 'Matsumoto et al. 2005, ApJ, 626, 1, pp. 31-43',
            'year': 2004,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [1.4300, 1.5300, 1.6300, 1.7300, 1.8300, 1.9300, 2.0300, 2.1400, 2.2400, 2.3400, 2.4400, 2.5400, 2.8800, 2.9800, 3.0700, 3.1700, 3.2800, 3.3800, 3.4800, 3.5800, 3.6800, 3.7800, 3.8800, 3.9800],
                'ebl': [70.1000, 71.3000, 65.9000, 58.7000, 51.0000, 43.2000, 39.2000, 35.4000, 29.7000, 24.4000, 22.7000, 23.7000, 19.5000, 18.6000, 19.7000, 16.4000, 13.6000, 14.6000, 14.4000, 15.5000, 13.1000, 13.5000, 15.3000, 15.5000],
                'ebl_err_low': [13.2000, 12.8000, 11.9000, 10.3000, 8.8000, 7.9000, 7.0000, 6.0000, 5.3000, 4.8000, 4.5000, 4.2000, 3.5000, 3.2000, 3.0000, 2.9000, 3.0000, 3.0000, 3.0000, 3.0000, 3.7000, 3.3000, 3.6000, 3.9000],
                'ebl_err_high': [13.2000, 12.8000, 11.9000, 10.3000, 8.8000, 7.9000, 7.0000, 6.0000, 5.3000, 4.8000, 4.5000, 4.2000, 3.5000, 3.2000, 3.0000, 2.9000, 3.0000, 3.0000, 3.0000, 3.0000, 3.7000, 3.3000, 3.6000, 3.9000]
                }
            }
        ,
        {
            'full_name': 'Matsuura et al. 2010 (AKARI)',
            'id': 'matsuura2010',
            'reference': 'ApJ, submitted',
            'year': 2010,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [65.0000, 90.0000, 140.0000, 160.0000],
                'ebl': [12.4600, 22.3300, 20.1400, 13.6900],
                'ebl_err_low': [9.2300, 4.6600, 3.4200, 3.9300],
                'ebl_err_high': [9.2300, 4.6600, 3.4200, 3.9300]
                }
            }
        ,
        {
            'full_name': 'Mattila 1990',
            'id': 'mattila1990',
            'reference': 'Hauser & Dwek 2001',
            'year': 1990,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [0.4000],
                'ebl': [46.0000],
                'ebl_err_low': [0.0000],
                'ebl_err_high': [0.0000]
                }
            }
        ,
        {
            'full_name': 'Metcalfe et al. 2003 (ISO)',
            'id': 'metcalfe2003',
            'reference': 'A&A, 407',
            'year': 2003,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [7.0000, 15.0000],
                'ebl': [0.4900, 2.7000],
                'ebl_err_low': [0.2000, 0.6200],
                'ebl_err_high': [0.2000, 0.6200]
                }
            }
        ,
        {
            'full_name': 'Papovich et al. 2004 (SPITZER)',
            'id': 'papovich2004',
            'reference': 'ApJS, 154, 1, pp. 70-74',
            'year': 2004,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [24.0000],
                'ebl': [2.7000],
                'ebl_err_low': [0.7000],
                'ebl_err_high': [1.1000]
                }
            }
        ,
        {
            'full_name': 'Thompson et al. 2007 (NICMOS)',
            'id': 'thompson2007',
            'reference': 'ApJ, 657, 2, 669-680',
            'year': 2007,
            'is_lower_limit': 1,
            'is_upper_limit': 0,
            'data': {
                'lambda': [1.1000, 1.6000],
                'ebl': [6.3000, 6.9000],
                'ebl_err_low': [0.3000, 0.3000],
                'ebl_err_high': [3.0000, 3.0000]
                }
            }
        ,
        {
            'full_name': 'Toller 1983/Leinert 1998',
            'id': 'toller1983',
            'reference': 'Hauser & Dwek 2001',
            'year': 1983,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [0.4400],
                'ebl': [60.0000],
                'ebl_err_low': [0.0000],
                'ebl_err_high': [0.0000]
                }
            }
        ,
        {
            'full_name': 'Wright & Reese 2000 (DIRBE)',
            'id': 'wrightReese2001',
            'reference': '',
            'year': 2000,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [1.2500, 2.2000, 3.5000],
                'ebl': [29.0000, 20.0000, 12.0000],
                'ebl_err_low': [16.0000, 6.0000, 3.0000],
                'ebl_err_high': [16.0000, 6.0000, 3.0000]
                }
            }
        ,
        {
            'full_name': 'Keenan et al. 2010',
            'id': 'keenan2010',
            'reference': 'The Astrophysical Journal, Volume 723, Issue 1, pp. 40-46 (2010)',
            'year': 2010,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [1.25, 1.6, 2.12],
                'ebl': [11.7000, 11.5000, 10.0000],
                'ebl_err_low': [2.6000, 1.5000, 0.8000],
                'ebl_err_high': [2.6000, 1.5000, 0.8000]
                }
            }
        ,
        {
            'full_name': 'Bethermin et al. 2010 (SPITZER)',
            'id': 'bethermin2010',
            'reference': 'Astronomy and Astrophysics, Volume 512, id.A78',
            'year': 2010,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [24., 70., 160.],
                'ebl': [2.86, 6.6, 14.6],
                'ebl_err_low': [0.16, 0.6, 2.9],
                'ebl_err_high': [0.19, 0.7, 7.1]
                }
            }
        ,
        {
            'full_name': 'Bethermin et al. 2010 (LL, SPITZER)',
            'id': 'bethermin2010ll',
            'reference': 'Astronomy and Astrophysics, Volume 512, id.A78',
            'year': 2010,
            'is_lower_limit': 1,
            'is_upper_limit': 0,
            'data': {
                'lambda': [24., 70., 160.],
                'ebl': [2.29, 5.4, 8.9],
                'ebl_err_low': [0.09, 0.4, 1.1],
                'ebl_err_high': [0.09, 0.4, 1.1]
                }
            }
        ,
        {
            'full_name': 'Berta et al. 2010 (Herschel/PEP)',
            'id': 'berta2010',
            'reference': 'Astronomy and Astrophysics, Volume 518, id.L30',
            'year': 2010,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [100., 160.],
                'ebl': [6.36, 6.58],
                'ebl_err_low': [1.67, 1.62],
                'ebl_err_high': [1.67, 1.62]
                }
            },
        {
            'full_name': 'Matsuoka et al. 2011 (Pioneer 10/11)',
            'id': 'matsuoka2011',
            'reference': 'The Astrophysical Journal, Volume 736, Issue 2, article id. 119',
            'year': 2011,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [.44, .64],
                'ebl': [7.9, 7.7],
                'ebl_err_low': [4., 5.8],
                'ebl_err_high': [4., 5.8]
                }
            },
        {
            'full_name': 'Mattila et al. 2011',
            'id': 'matilla2011',
            'reference': 'Invited talk, IAU Symposium No.284',
            'year': 2011,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [.4],
                'ebl': [7.2],
                'ebl_err_low': [2.],
                'ebl_err_high': [4.]
                }
            },
        {
            'full_name': 'Mattila et al. 2011 (UL)',
            'id': 'matilla2011ul',
            'reference': 'Invited talk, IAU Symposium No.284',
            'year': 2011,
            'is_lower_limit': 0,
            'is_upper_limit': 1,
            'data': {
                'lambda': [.52],
                'ebl': [12.],
                'ebl_err_low': [0.],
                'ebl_err_high': [0.]
                }
            },
        {
            # EBL(160m) = 0.77+/-0.04+/-0.12 MJy/sr
            # EBL(100m) = 0.24+/-0.08+/-0.04 MJy/sr
            'full_name': 'Penin et al. 2012',
            'id': 'penin2012',
            'reference': 'Astronomy & Astrophysics, Volume 543, id.A123',
            'year': 2012,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [100., 160.],
                'ebl': [3000. / 100. * .24, 3000. / 160. * .77],
                'ebl_err_low': [3000. / 100. * .08, 3000. / 160. * .12],
                'ebl_err_high': [3000. / 100. * .08, 3000. / 160. * .12]
                }
            },
        {
            'full_name': 'Xu et al. 2005 (GALEX)',
            'id': 'xu2005',
            'reference': 'The Astrophysical Journal, Volume 619, Issue 1, pp. L11-L14.',
            'year': 2005,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [.153, .2310],
                'ebl': [1.03, 2.25], #[.68, .99],
                'ebl_err_low': [.15, .32],#[.1, .15],
                'ebl_err_high':[.15, .32]# [.1, .15]
                }
            },
         {
            'full_name': 'Gardner et al. 2000 (STIS)',
            'id': 'gardner2000',
            'reference': 'The Astrophysical Journal, Volume 542, Issue 2, pp. L79-L82.',
            'year': 2000,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [.1595, .2365],
                'ebl': [2.9, 3.6], #[.68, .99],
                'ebl_err_low': [.4, .5],#[.1, .15],
                'ebl_err_high':[.6, .7]# [.1, .15]
                }
            },
         {
            'full_name': 'Voyer et al. 2011 (SBC/GALEX)',
            'id': 'voyer2011',
            'reference': 'The Astrophysical Journal, Volume 736, Issue 2, article id. 80 (2011).',
            'year': 2011,
            'is_lower_limit': 0,
            'is_upper_limit': 0,
            'data': {
                'lambda': [.1572, .1572],
                'ebl': [1.3, 1.6],
                'ebl_err_low': [.2, .2],
                'ebl_err_high':[.2, .2]
                }
            }
       ]
    return ebl_limit_collection

#---------------------------------------------------------------------------
def plot_ebl_measurement_collection(color='.75', cm=None, yearmin=None, yearmax=None, nolabel=False) :

    ebl_m = get_ebl_measurement_collection()
    ebl_m = sorted(ebl_m, key=lambda t: t['year'])

    if yearmin or yearmax:
        ebl_m_n = []
        for m in ebl_m :
            if yearmin and m['year'] >= yearmin :
                ebl_m_n.append(m)
            if yearmax and m['year'] < yearmax :
                ebl_m_n.append(m)
        ebl_m = ebl_m_n

    t = 0.
    markers = ['*', '<', '>', 'H', '^', 'd', 'h', 'o', 'p', 's', 'v']

    for mi, m in enumerate(ebl_m):
        yerr = [m['data']['ebl_err_low'], m['data']['ebl_err_high']]
        if m['is_upper_limit'] or m['is_lower_limit'] :
            yerr = [ 10. ** (np.log10(x) + 0.14) - x for x in  m['data']['ebl']]
            if m['is_upper_limit'] :
                yerr = [yerr, map(lambda x: 0., range(len(yerr)))]
            else :
                yerr = [map(lambda x: 0., range(len(yerr))), yerr]                
        if cm :
            color = cm(float(mi) / (len(ebl_m) - 1.))
        label = m['full_name']
        if nolabel : label=""
        plt.errorbar(x=m['data']['lambda'], y=m['data']['ebl'],
                     yerr=yerr, label=label,
                     marker=markers[mi % len(markers)], color=color, mec=color,
                     linestyle='None', lolims=m['is_upper_limit'], uplims=m['is_lower_limit'])

#===========================================================================

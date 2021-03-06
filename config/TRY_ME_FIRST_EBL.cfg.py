{
    'plottitle' : r'Test PopII',
    'models' :
    [
        {
            'name' : 'PopII fesc=0.',
            # BPL 10. ** (z+1)
            'sfrf' : 'lambda p, x : p[0] * (10. ** (x - p[2])) ** -p[1] * ((1. + (10. ** (x - p[2])) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            'sfrp0' : [1E-1, -1., 2.1, .1, .6],
            'spstype' : 'SB99',
            'spsfile' : 'data/ssp/starburst99/sb99default.spectrum1.gz',# Starburst 99 test
            'lmin' : 1E-2,
            'lmax' : 3E3,
            'lfsteps' : 201,#201
            'zmin' : 0.,
            'zmax' : 10.,
            'zintmax' : 35.,
            'zsteps' : 201,#201
            'ltintsteps' : 201,#201
            'eblzmin' : 0.,
            'eblzmax' : 10.,
            'eblzsteps' : 101,#
            'cosmo' : [.7, .3, .7],
            'ebv' : .15,
            'ir_wave_start' : 5.5,#5.5
            'ir_fac' : 3E9,
            'dust' : True,
            'nebula' : True,
            'fesc' : 0.            
            }
        ]
    }

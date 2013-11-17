{
    'plottitle' : r'PopIII test',
    'models' :
    [
        {
            'name' : 'PopIII fesc=.5',
            # BPL 10. ** (z+1)
            'sfrf' : 'lambda p, x : p[0] * (10. ** (x - p[2])) ** -p[1] * ((1. + (10. ** (x - p[2])) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            'sfrp0' : [2E-4, -2., 7.5, 0.15, 0.2],
            'spsfile' : '/Users/mraue/Stuff/unix/python/a3p2/data/ssp/tumlinson2006/ssp_Z=0_caseA_v2.dat',
            'spstype' : 'TUM2006',
            'lmin' : 0.0001,
            'lmax' : 1E3,
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
            'dust' : False,
            'nebula' : True,
            'fesc' : 0.
            }, 
        {
            'name' : 'PopIII fesc=0.',
            # BPL 10. ** (z+1)
            'sfrf' : 'lambda p, x : p[0] * (10. ** (x - p[2])) ** -p[1] * ((1. + (10. ** (x - p[2])) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            'sfrp0' : [2E-4, -2., 7.5, 0.15, 0.2],
            'spsfile' : '/Users/mraue/Stuff/unix/python/a3p2/data/ssp/tumlinson2006/ssp_Z=0_caseA_v2.dat',
            'spstype' : 'TUM2006',
            'lmin' : 0.0001,
            'lmax' : 1E3,
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
            'dust' : False,
            'nebula' : True,
            'fesc' : .5
            },
       {
            'name' : 'PopIII fesc=1.',
            # BPL 10. ** (z+1)
            'sfrf' : 'lambda p, x : p[0] * (10. ** (x - p[2])) ** -p[1] * ((1. + (10. ** (x - p[2])) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            'sfrp0' : [2E-4, -2., 7.5, 0.15, 0.2],
            'spsfile' : '/Users/mraue/Stuff/unix/python/a3p2/data/ssp/tumlinson2006/ssp_Z=0_caseA_v2.dat',
            'spstype' : 'TUM2006',
            'lmin' : 0.0001,
            'lmax' : 1E3,
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
            'dust' : False,
            'nebula' : True,
            'fesc' : 1.
            }
        ]
    }

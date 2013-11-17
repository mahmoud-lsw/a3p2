{
    'plottitle' : r'Best fit - T2007',
    'models' :
    [
        {
            'name' : 'PopII',
            #'sfrf' : 'lambda p, x : np.where(x < p[0], p[1] * (x / p[0]) ** -p[2],  p[1] * (x / p[0]) ** -p[3])',
            #'sfrp0' : [2.2, .07, -3.5, 1.2],
            # Cole et al. 2001
            #'sfrf' : 'lambda p, x : (p[0] + (x - 1.) * p[1]) / (1. + ((x - 1.) / p[2]) ** p[3]) * .7',
            #'sfrp0' : [0.0166, 0.1848, 1.9474, 2.6316],
            #'sfrp0' : [0.0, 0.0798, 1.658, 3.105],
            #'sfrp0' : [0.0166, 0.3, 2., 4.],
            #'sfrp0' : [0.017, 0.13, 3.3, 5.3],
            # BPL
            #'sfrf' : 'lambda p, x : p[0] * (x / p[2]) ** -p[1] * ((1. + (x / p[2]) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            #'sfrp0' : [.16, -3.5, 2.2, 10., 3.5],
            # BPL 10. ** (z+1)
            'sfrf' : 'lambda p, x : p[0] * (10. ** (x - p[2])) ** -p[1] * ((1. + (10. ** (x - p[2])) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            'sfrp0' : [1E-1, -1.5, 2.1, .3, .6],
            'spsfile' : '/Users/mraue/Stuff/work/ebl/eblmodel-ir-new/bruzual_charlot_2003/bc03/models/Padova1994/salpeter/bc2003_lr_m62_salp_ssp.ised_ASCII.gz',# BC 2003, Salpeter IMF, Z=Z_solar
            'lmin' : 0.08,
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
            'nebula' : True
            },
        {
            'name' : 'PopIILM',
            #'sfrf' : 'lambda p, x : np.where(x < p[0], p[1] * (x / p[0]) ** -p[2],  p[1] * (x / p[0]) ** -p[3])',
            #'sfrp0' : [2.2, .07, -3.5, 1.2],
            # Cole et al. 2001
            #'sfrf' : 'lambda p, x : (p[0] + (x - 1.) * p[1]) / (1. + ((x - 1.) / p[2]) ** p[3]) * .7',
            #'sfrp0' : [0.0166, 0.1848, 1.9474, 2.6316],
            #'sfrp0' : [0., 0.0798, 1.658, 3.105],
            #'sfrp0' : [1E-10, 0.01, 6.3, 5.3],
            # BPL
            #'sfrf' : 'lambda p, x : p[0] * (x / p[2]) ** -p[1] * ((1. + (x / p[2]) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            #'sfrp0' : [.14, -4., 4., 5., 10.],            
            # BPL 10. ** (z+1)
            'sfrf' : 'lambda p, x : p[0] * (10. ** (x - p[2])) ** -p[1] * ((1. + (10. ** (x - p[2])) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            'sfrp0' : [7E-3, -2., 3., .25, .4],
            # Parabola 10. ** x
            #'sfrf' : 'lambda p, x : 10. ** (p[0] + p[1] * (x - p[2]) + p[3] * (x - p[2]) ** 2.)',
            #'sfrp0' : [-1., 2., 5., -1.],
            'spsfile' :'/Users/mraue/Stuff/work/ebl/eblmodel-ir-new/bruzual_charlot_2003/bc03/models/Padova1994/salpeter/bc2003_lr_m42_salp_ssp.ised_ASCII',# BC 2003 salpeter, z=0.004
            'lmin' : 0.08,
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
            'nebula' : True
            },
        {
            'name' : 'PopIII',
            # BPL 10. ** (z+1)
            'sfrf' : 'lambda p, x : p[0] * (10. ** (x - p[2])) ** -p[1] * ((1. + (10. ** (x - p[2])) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            'sfrp0' : [.9E-5, -.6, 7.5, .3, .25],
            'spsfile' : '/Users/mraue/Stuff/unix/python/a3p2/data/ssp/tumlinson2006/ssp_Z=0_caseA_v2.dat',
            'spstype' : 'TUM2006',
            'lmin' : 0.08,
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
            'dust' : False,
            'nebula' : True
            }
        ]
    }

{
    'plottitle' : r'Paper 1 - v1',
    'scansfrd' : {
        'beta_a' : [.3, .7, 1.],
        'zpeakmin' : .2,
        'zpeakmax' : 2.4,
        'zpeakn' : 14,#14
        'lrho0min' : -1.2,
        'lrho0max' : -.5,
        'lrho0n' : 8,#8
        'rhoz0' : .012        
        },
    'models' :
    [
        {
            'name' : 'Pop2',
            'sfrf' : 'lambda p, x : np.where(x < p[4], np.where(x < p[0], p[1] * (x / p[0]) ** -p[2],  p[1] * (x / p[0]) ** -p[3]), 0.)',
            'sfrp0' : [2.2, .09, -3., .3, 5.],
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
            #'sfrf' : 'lambda p, x : p[0] * (10. ** (x - p[2])) ** -p[1] * ((1. + (10. ** (x - p[2])) ** p[3]) / 2.) ** ((p[1] - p[4]) / p[3])',
            #'sfrp0' : [1E-1, -1., 2.1, .1, .6],
            'spsfile' : '/Users/mraue/Stuff/work/ebl/eblmodel-ir-new/bruzual_charlot_2003/bc03/models/Padova1994/salpeter/bc2003_lr_m62_salp_ssp.ised_ASCII.gz',# BC 2003, Salpeter IMF, Z=Z_solar
            'lmin' : 0.01,
            'lmax' : 3E3,
            'lfsteps' : 201,#201
            'zmin' : 0.,
            'zmax' : 10.,
            'zintmax' : 10.,
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
            'fesc' : .0
            }
        ]
    }

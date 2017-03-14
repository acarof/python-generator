#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy as np
import imp

from utils import *


def main(inputs, paths):
    print """
    """

    task = {
        'KIND_RUN' : 'TEST_JACOB_RABI',
         'TEST' : 'YES',
        'FILE_INIT': 'initial',
        'MOL_NAME': 'ETHYLENE',
        'NATOM_MOL': 6,
        'NATOMS'   : 12,
        'VECTA': [3.527, 0.784, -0.166],
        'VECTB': [0.0, 0.0, 0.0],
        'VECTC': [0.0, 0.0, 0.0],
        'SIZE_CRYSTAL': [2, 1, 1],
        'COORD_CHARGE': [2, 1, 1],
        'FIRST_DIABAT': 2,
        'RCUT': 12,
        'SYSTEM': 'CRYSTAL'
    }
    inputs.update(task)



    # SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
    bucket = Bucket(inputs)
    bucket.name()
    paths.update({'bucket': bucket.path})

    task = Dir(inputs.get('INPUT_INFO'))
    paths.update( {'task' : task.path} )

    templates = Dir('templates', paths)
    templates.checkdir()
    templates.clean()

    bin = Dir('bin', paths)
    bin.checkdir()

    # FIND CP2K PATHS
    try:
        local_paths = Dir('local_paths', paths)
        local_paths.checkdir()
        cp2k_file = open(paths.get('local_paths') + 'cp2k.path', 'r')
        paths.update({'cp2k': cp2k_file.read().rstrip()})
        if not os.path.isfile(paths.get('cp2k')):
            raise SystemExit('WARNING: check path for CP2K executable in local_paths/cp2k.path')
    except:
        raise SystemExit("WARINING: please provide the path for CP2K executable in local_paths/cp2k.path")

    system = inputs.get('SYSTEM')
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    print paths.get('bucket') + '/scripts/utils_scripts.py'
    scripts = imp.load_source('scripts', paths.get('bucket') + '/scripts/utils_scripts.py')


    # 1. CREATE INITIAL DIRECTORY

    initial = Dir('initial', paths)
    initial.rm_mkdir()
    for i in range(100):
        cmd = 'cp %sCOORD/COORD_%d initial/pos-%d.init' % ( paths.get('task'), i+1, i )
        os.system(cmd)
        file = open('initial/vel-%d.init' % i, 'w')
        results = ''
        for atom in range(inputs.get('NATOMS')):
            results += '0.00  0.00  0.00\n'
        file.write(results)
        file.close()
    ndir = 0



    # 2. CHECK HAMILTONIAN, NACV and NACE

    subtask = {
        'TEMPLATE_FILE': 'FSSH_CORE.template',
        'FORCEFIELD_FILE': 'FSSH_FF.template',
        'NUMBER_INIT' : 100,
        'STEPS' : 3,
        'PRINT' : 1,
        'SCALING' : 0.065190,
        'PROPAGATION' : 'FSSH',
        'TIMESTEP'    : 0.01
    }
    final_step = subtask.get('STEPS')*subtask.get('TIMESTEP')
    couplings = []
    nace = []
    delta_e = []
    nacv = []
    for init in range(subtask.get('NUMBER_INIT')):
        config = Config(inputs, paths, INIT=init, **subtask)
        print "GO FOR RUN %d" % ndir
        dir = scripts.FSSHRun( 'run-%d' % ndir)
        ndir = config.run(ndir)
        if os.path.exists('run-%d' % (ndir - 1)):
            pass
            #os.system('rm -rf run-%d' % (ndir - 1))
        else:
            print " ERROR IN CP2K FOR THOSE PARAMETERS:"
            print dict
            print " TO BE COMPARED WITH:"
            print dict_prev
            sys.exit()
        dict_prev = dict
        os.chdir('run-%d' % (ndir -1))
        couplings.append( dir.extract('Couplings').get( final_step )[0] * 27.211399 ) # Coupling in ev (from Ha)
        nace.append( dir.extract('NACE').get( final_step)[0] )
        for atom in range(1,1+inputs.get('NATOMS')):
            for dim in range(1,4):
                nacv.append( dir.extract('NACV').get(final_step).get(atom).get(dim)[0][1] )
        for atom in range(1,1+inputs.get('NATOMS')):
            for dim in range(1,4):
                nacv.append( dir.extract('NACV').get(final_step).get(atom).get(dim)[1][0] )
        delta_e.append( dir.extract('Delta_E').get( final_step)[0] * 27.211399  )
        os.chdir('..')



    file = open( paths.get('task') + 'DATA_COMPARE/' + 'Hab-laura-160419.dat')
    laura_coupling = [ float(line) * 0.04336412 for line in file.readlines()] # Coupling in eV (from kcal/mol)
    file1 = open( paths.get('task') + 'DATA_COMPARE/' + 'E1-laura-160419.dat')
    file2 = open( paths.get('task') + 'DATA_COMPARE/' + 'E2-laura-160419.dat')
    laura_delta_e = [ ( - float(line1.split()[1]) + float(line2.split()[1]) )*0.04336412
                      for line1, line2 in zip(file1.readlines(), file2.readlines())
                      ]
    file = open( paths.get('task') + 'DATA_COMPARE/' + 'nace-laura-160708.dat')
    laura_nace = [float(line)  for line in file.readlines()]
    file = open( paths.get('task') + 'DATA_COMPARE/' + 'nacv-laura-160708.dat')
    laura_nacv = [float(line)  for line in file.readlines()]


    print scripts.rmse( delta_e, laura_delta_e)
    print scripts.rmse( couplings, laura_coupling )
    print scripts.rmse( nace, laura_nace)
    print scripts. rmse( nacv, laura_nacv)




    # CHECK RABI OSCILLATION

    subtask = {
        'TEMPLATE_FILE': 'FSSH_CORE.template',
        'FORCEFIELD_FILE': 'FSSH_FF.template',
        'NUMBER_INIT': 3,
        'STEPS': 3,
        'PRINT': 1,
        'SCALING': 0.0065190,
        'PROPAGATION': 'FROZEN_HAMILTONIAN'
    }
    for init in range(subtask.get('NUMBER_INIT')):
        populations = []
        config = Config(inputs, paths, INIT=init, **subtask)
        print "GO FOR RUN %d" % ndir
        dir = scripts.FSSHRun( 'run-%d' % ndir)
        ndir = config.run(ndir)
        if os.path.exists('run-%d' % (ndir - 1)):
            pass
            #os.system('rm -rf run-%d' % (ndir - 1))
        else:
            print " ERROR IN CP2K FOR THOSE PARAMETERS:"
            print dict
            print " TO BE COMPARED WITH:"
            print dict_prev
            sys.exit()
        dict_prev = dict
        os.chdir('run-%d' % (ndir -1))

        pop = dir.extract('Populations')
        for time in sorted(pop):
            populations.append( pop.get(time)[1] )
        populations = populations[4:]
        coupling = dir.extract('Couplings').get( 0 )[0]
        delta_e = dir.extract('Delta_E').get( 0 )[0]

        rabi = scripts.rabi_oscillation( coupling, delta_e, 0.1, subtask.get('STEPS')*5  )

        print scripts.rmse(populations, rabi)
        os.chdir('..')

#    # CHECK DIAGONAL FORCES
#
#    subtask = {
#        'NUMBER_INIT': 100,
#        'STEPS': 3,
#        'PRINT': 1,
#        'SCALING': 0.0
#    }
#
#    dict_list = [
#        {
#            'TEMPLATE_FILE': 'FSSH_CORE.template',
#            'FORCEFIELD_FILE': 'FSSH_FF.template',
#            'PROPAGATION': 'BORN_OPPENHEIMER'
#        },
#        {
#            'TEMPLATE_FILE': 'FIST_TEMPLATE'
#        }
#    ]
#
#    for dict_ in dict_list:
#        for init in range(subtask.get('NUMBER_INIT')):
#            forces = []
#            subsubtask = subtask.update(dict_)
#            config = Config(inputs, paths, INIT=init, **subsubtask)
#            print "GO FOR RUN %d" % ndir
#            dir = scripts.FSSHRun('run-%d' % ndir)
#            ndir = config.run(ndir)
#            if os.path.exists('run-%d' % (ndir - 1)):
#                pass
#                # os.system('rm -rf run-%d' % (ndir - 1))
#            else:
#                print " ERROR IN CP2K FOR THOSE PARAMETERS:"
#                print dict
#                print " TO BE COMPARED WITH:"
#                print dict_prev
#                sys.exit()
#            dict_prev = dict
#            os.chdir('run-%d' % (ndir - 1))
#
#            os.chdir('..')
#


    # CHECK FORCES ANALYTICS

    subtask = {
        'TEMPLATE_FILE': 'FSSH_CORE.template',
        'FORCEFIELD_FILE': 'FSSH_FF.template',
        'NUMBER_INIT': 1,
        'STEPS': 200,
        'PRINT': 1,
        'SCALING': 0.0065190,
        'PROPAGATION': 'BORN_OPPENHEIMER',
        'ANALYTICS' : 'T'
        }

    list_state = [1, 2]
    forces = []
    exforces = []
    for state in list_state:
        for init in range(subtask.get('NUMBER_INIT')):
                   subtask.update({'FIRST_DIABAT' : state})
                   config = Config(inputs, paths, INIT=init, **subtask)
                   print "GO FOR RUN %d" % ndir
                   dir = scripts.FSSHRun('run-%d' % ndir)
                   ndir = config.run(ndir)
                   if os.path.exists('run-%d' % (ndir - 1)):
                       pass
                       # os.system('rm -rf run-%d' % (ndir - 1))
                   else:
                       print " ERROR IN CP2K FOR THOSE PARAMETERS:"
                       print dict
                       print " TO BE COMPARED WITH:"
                       print dict_prev
                       sys.exit()
                   dict_prev = dict
                   os.chdir('run-%d' % (ndir - 1))

                   frc = dir.extract('Forces')
                   for time in sorted(frc):
                       for i in range(len(frc.get(time))):
                           if i%3 != 0:
                               for element in frc.get(time)[i]:
                                    forces.append( element)

                   frc = dir.extract('Forces', filename = 'run-exfrc-1.xyz')
                   for time in sorted(frc):
                       for i in range(len(frc.get(time))):
                           if i%3 != 0:
                               for element in frc.get(time)[i]:
                                    forces.append( element)


                   os.chdir('..')

    print scripts.rmse(forces, exforces)






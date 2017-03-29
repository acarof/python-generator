#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy as np
import imp
import subprocess
from scipy.optimize import curve_fit



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


    system = inputs.get('SYSTEM')
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
        from utils import CP2KOSFIST as ConfigFIST
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
            results += '0.00 0.00  0.00  0.00\n'
        file.write(results)
        file.close()
    ndir = 0
    final = "\n\n\n\n"
    final += " RESULTS OF TESTS JACOB_RABI\n\n"
    final += " Date: %s \n\n" % time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    final += "CP2K Version:\n"
    version = subprocess.check_output([ paths.get('cp2k'), '--version'])
    final += version
    final += "\n\n"



    # 2. CHECK HAMILTONIAN, NACV and NACE

    subtask = {
        'TEMPLATE_FILE': 'FSSH_CORE.template',
        'FORCEFIELD_FILE': 'FSSH_FF.template',
        'NUMBER_INIT' : 5 ,
        'STEPS' : 3,
        'PRINT' : 1,
        'SCALING' : 0.065190,
        'PROPAGATION' : 'FSSH',
        'TIMESTEP'    : 0.01
    }
    final_step = subtask.get('STEPS')*subtask.get('TIMESTEP')
    couplings = []
    #nace = []
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
        #nace.append( dir.extract('NACE').get( final_step)[0] )
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
    #file = open( paths.get('task') + 'DATA_COMPARE/' + 'nace-laura-160708.dat')
    #laura_nace = [float(line)  for line in file.readlines()]
    file = open( paths.get('task') + 'DATA_COMPARE/' + 'nacv-laura-160708.dat')
    laura_nacv = [float(line)  for line in file.readlines()]


    final += "1. CHECK HAMITLONIAN, NACE AND NACV\n"

    final += "Number of config: %s\n " % subtask.get('NUMBER_INIT')
    final += "Number of steps:  %s\n\n" % subtask.get('STEPS')


    final += "RMSE between CP2K couplings and NAMD-FSSH couplings " \
          " is: %s eV \n" \
          % (scripts.rmse( couplings, laura_coupling ))
    final += "RMSE between CP2K Delta_E and NAMD-FSSH Delta_E " \
          " is: %s eV \n" \
          % (scripts.rmse( delta_e, laura_delta_e))
    #final += "RMSE between CP2K NACE and NAMD-FSSH NACE " \
    #      " is: %s atomic unit \n" \
    #      % (scripts.rmse( nace, laura_nace ))
    final += "RMSE between CP2K NACV and NAMD-FSSH NACV " \
          " is: %s atomic unit \n" \
          % (scripts.rmse( nacv, laura_nacv ))
    final += "\n\n"


    # CHECK DECO DECAY: FROZEN HAMILTONIN SCALING=0
    
    subtask = {
        'TEMPLATE_FILE': 'FSSH_CORE_without_constraint.template',
        'FORCEFIELD_FILE': 'FSSH_FF_without_constraint.template',
        'NUMBER_INIT': 2,
        'STEPS': 50,
        'PRINT': 1,
        'SCALING': 0.0,
        'PROPAGATION': 'FROZEN_HAMILTONIAN',
        'DECO': 'DAMPING',
        'EDC_C': 1.0,
        'EDC_E0': 0.0
    }
     

    rmse = []
    tau_calculated = []
    tau_fitted = []
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

        step = []
        pop = dir.extract('Populations')
        for times in sorted(pop):
            populations.append( pop.get(times)[1] )
            step.append(times)
        populations = np.array(populations[4:])
        if populations[0] < 0.2: 
            populations = 1 - populations
        else:
            pass

        step = np.array(step[4:])
        def func(x1, a, b):
            return a * np.exp(-2*b*x1)
        popt, pcov = curve_fit(func, step, populations)
        tau_fit = 1/popt[1]
        
         
        delta_e = dir.extract('Delta_E').get( 0 )[0]
        au_to_fs = 41.34137
        tau_cal = 1/(delta_e*au_to_fs)
        diff = abs(tau_cal- tau_fit)
        tau_calculated.append(tau_cal)
        x = np.array(tau_calculated)
        tau_fitted.append(tau_fit)
        y = np.array(tau_fitted)


        rmse.append(scripts.rmse(x, y))
        os.chdir('..')
    final += "2. CHECK DECO DECAY: FROZEN HAMILTONIN SCALING=0 \n"
    
    final += "Number of config: %s\n " % subtask.get('NUMBER_INIT')
    final += "Number of steps:  %s\n\n" % subtask.get('STEPS')
   # final += "difference in decay time is: %s\n" % diff
    rmse_mean = np.mean(rmse)
    final += "RMSE between calculated decoherence time and fitted decoherence time" \
          " is: %s\n" \
          % rmse_mean
    final += "\n\n"



    # CHECK RABI OSCILLATION

    subtask = {
        'TEMPLATE_FILE': 'FSSH_CORE_without_constraint.template',
        'FORCEFIELD_FILE': 'FSSH_FF_without_constraint.template',
        'NUMBER_INIT': 1,
        'STEPS': 5,
        'PRINT': 1,
        'SCALING': 0.0065190,
        'PROPAGATION': 'FROZEN_HAMILTONIAN',
        'DECO' : 'NO_DECO_CORR',
        'TIMESTEP' : 0.5
    }
    rmse = []
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
        for times in sorted(pop):
            populations.append( pop.get(times)[1] )
        populations = populations[4:]
        coupling = dir.extract('Couplings').get( 0 )[0]
        delta_e = dir.extract('Delta_E').get( 0 )[0]

        rabi = scripts.rabi_oscillation(  coupling, delta_e, subtask.get('TIMESTEP') / 5, subtask.get('STEPS')*5  )

        rmse.append(scripts.rmse(populations, rabi))
        os.chdir('..')


    final += "3. CHECK RABI OSCILLATION\n"

    final += "Number of config: %s\n " % subtask.get('NUMBER_INIT')
    final += "Number of steps:  %s\n\n" % subtask.get('STEPS')
    rmse_mean = np.mean( rmse)
    final += "RMSE between CP2K Population and Rabi Oscillation " \
          " is: %s\n" \
          % rmse_mean
    final += "\n\n"


    # CHECK DIAGONAL FORCES

    subtask = {
        'NUMBER_INIT': 100,
        'STEPS': 2,
        'PRINT': 1,
        'SCALING': 0.0
    }

    dict_list = [
        {
            'TEMPLATE_FILE': 'FSSH_CORE.template',
            'FORCEFIELD_FILE': 'FSSH_FF.template',
            'PROPAGATION': 'BORN_OPPENHEIMER'
        },
        {
            'TEMPLATE_FILE': 'FIST_TEMPLATE',
            'PROPAGATION'  : 'FIST'
        }
    ]

    forces = [ [], []]
    ind = 0
    for dict_ in dict_list:
        for init in range(subtask.get('NUMBER_INIT')):
            subsubtask = subtask
            subsubtask.update(dict_)
            if subsubtask.get('PROPAGATION') == 'FIST':
                config = ConfigFIST(inputs, paths, INIT=init, **subsubtask)
            else:
                config = Config(inputs, paths, INIT=init, **subsubtask)
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
            for times in sorted(frc):
                for i in range(len(frc.get(times))):
                    for element in frc.get(times)[i]:
                        forces[ind].append(element)
            os.chdir('..')
        ind += 1

    final += "4. CHECK DIAGONAL FORCES\n"

    final += "Number of config: %s\n " % subtask.get('NUMBER_INIT')
    final += "Number of steps:  %s\n\n" % subtask.get('STEPS')
    final += "RMSE between diagonal forces (BO) and FIST forces " \
          " is: %s atomic unit \n" \
          % (scripts.rmse( forces[0], forces[1]) )
    final += "\n\n"


    # CHECK FORCES ANALYTICS

    subtask = {
        'TEMPLATE_FILE': 'FSSH_CORE.template',
        'FORCEFIELD_FILE': 'FSSH_FF.template',
        'NUMBER_INIT': 1,
        'STEPS': 200,
        'PRINT': 1,
        'SCALING': 0.0065190,
        'PROPAGATION': 'BORN_OPPENHEIMER',
        'ANALYTICS' : 'T',
        'CC_CHARGED' : 1.369
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
                   for times in sorted(frc):
                       for i in range(len(frc.get(times))):
                           if i%3 != 0:
                               for element in frc.get(times)[i]:
                                    forces.append( element)

                   frc = dir.extract('Forces', filename = 'run-exfrc-1.xyz')
                   for times in sorted(frc):
                       for i in range(len(frc.get(times))):
                           if i%3 != 0:
                               for element in frc.get(times)[i]:
                                    forces.append( element)


                   os.chdir('..')

    final += "5. CHECK FORCES FORMULA\n"

    final += "Number of config: %s\n " % subtask.get('NUMBER_INIT')
    final += "Number of steps:  %s\n\n" % subtask.get('STEPS')
    final += "RMSE between FSSH forces (BO) and exact forces " \
          " is: %s atomic unit \n" \
          % scripts.rmse(forces, exforces)
    final += "\n\n"

    print final





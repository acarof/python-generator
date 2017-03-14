#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

def main(inputs, paths):
    print """
    """

    task = {
        'KIND_RUN' : 'TEST_CP2K',
        'TEMPLATE_FILE' : 'FSSH_CORE.template',
        'FORCEFIELD_FILE' : 'FSSH_FF.template',
         'TEST' : 'YES'
    }
    inputs.update(task)


    list_propagation = ['FSSH','BORN_OPPENHEIMER', 'TEST_HOP','FROZEN_HAMILTONIAN','CLASSICAL_PATH','GALILEAN']
    #list_propagation = ['FSSH']
    list_decoherences = ['NO_DECO_CORR', 'INSTANT_COLLAPSE', 'DAMPING', 'TRESH_ONLY_COLLAPSE']

    list_analytics   = ['T', 'F']
    #list_analytics = ['T']

#    list_collapse = ['T', 'F']
    #list_analytics = ['T']
    list_com = ['T', 'F']
    #list_collapse = ['T']
    list_first_diabat = [1, 2]
    #list_first_diabat = [1]
    list_rescaling    = ['SIMPLE', 'NACV']
    #list_rescaling = ['NACV']
    list_nacv         = ['TEST', 'CONTRIBUTION','TOTAL','FAST']
    #list_nacv         = ['TOTAL']
    list_reversal     = ['NEVER','ALWAYS','TRHULAR','SUBOTNIK']
    #list_reversal = ['ALWAYS']

    mega_list = [ { 'PROPAGATION' : prop,
                   # 'COLLAPSE'    : collapse,
                    'DECOHERENCES': deco,
                    'ANALYTICS'   : analytics,
                    'FIRST_DIABAT': diabat,
                    'METHOD_RESCALING' : rescaling,
                    'METHOD_ADIAB_NACV' : nacv,
                    'METHOD_REVERSAL'   : reversal,
                    'CENTER_OF_MASS'    : com
                   }
                  for prop in list_propagation
                 # for collapse in list_collapse
                  for deco in list_decoherences
                  for analytics in list_analytics
                  for diabat    in list_first_diabat
                  for rescaling in list_rescaling
                  for nacv      in list_nacv
                  for reversal  in list_reversal
                  for com       in list_com
                  ]

    systems = ['dimer', 'trimer', 'dimer_solvent']



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


    ndir= 0

    for system in systems:
        system_input = InputFile( paths.get('task') + system + '/input')
        os.system(' cp -r %s/initial/ %s' % (paths.get('task') + system, paths.get('bucket')))
        initial = Dir( 'initial', paths)
        initial.checkdir()     
        paths.update({'initial' :  initial.path})
        inputs.update(system_input.dict)
        inputs.update({'STEPS': 2})
        dict_prev = {}
        for dict in mega_list:
            if dict.get('PROPAGATION') == 'TEST_HOP':
               inputs.update({'STEPS':1})
            if system_input.dict.get('FILE_INIT') != 'initial_dimer':
               dict.update({'ANALYTICS' : 'F' })

            system = system_input.dict.get('SYSTEM')
            if system == 'CRYSTAL':
                from utils import CP2KOSFSSH as Config
            elif system == 'SOLVENT':
                from utils import CP2KOSwSolventFSSH as Config
            else:
                sys.exit()

            config = Config( inputs, paths, INIT = 1, **dict)
            print "GO FOR RUN %d" % ndir
            ndir = config.run(ndir)
            if os.path.exists('run-%d' % (ndir -1) ):
            #    pass
                os.system('rm -rf run-%d' % (ndir - 1))
            else:
                print " ERROR IN CP2K FOR THOSE PARAMETERS:"
                print dict
                print " TO BE COMPARED WITH:"
                print dict_prev
                sys.exit()
            dict_prev = dict




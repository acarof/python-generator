#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *
from multiprocessing import Pool




def run_fssh( dict_):
    inputs = dict_.get('INPUTS_DICT')
    paths  = dict_.get('PATHS_DICT')
    system = dict_.get('SYSTEM')
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    config = Config(inputs, paths, **dict_)
    ndir = dict_['NDIR']
    print "GO FOR RUN %d" % ndir
    config.run(ndir)


def main(inputs, paths):
    print """
    """


    task = {
        'KIND_RUN' : 'TASK234-REVERSAL',
        'STEPS'    : 20,
        'SELECT'   : 'T',
        'FIRST_ADIABAT' : 1,
        'NUMBER_CONFIG'        : 10,
        'NUMBER_REPEAT'  :  1
    }
    inputs.update(task)



    list_propagation = ['FSSH']
    list_decoherences = ['DAMPING']
    list_rescaling    = ['NACV']
    list_nacv         = ['FAST']
    list_reversal = ['NEVER', 'ALWAYS', 'TRHULAR', 'SUBOTNIK']
    #list_reversal = ['NEVER']
    list_init     = range(1, inputs.get('NUMBER_CONFIG') + 1)
    list_repeat   = range(inputs.get('NUMBER_REPEAT'))
    list_scaling = [0.0001, 0.001, 0.01, 0.03, 0.05, 0.07, 0.1]


    mega_list = [ { 'PROPAGATION' : prop,
                    'DECOHERENCES': deco,
                    'METHOD_RESCALING' : rescaling,
                    'METHOD_ADIAB_NACV' : nacv,
                    'METHOD_REVERSAL'   : reversal,
                    'INIT'              : init ,
                    'REPEAT'             : repeat,
                    'SCALING'            : scaling
                   }
                  for prop in list_propagation
                  for deco in list_decoherences
                  for rescaling in list_rescaling
                  for nacv      in list_nacv
                  for reversal  in list_reversal
                  for init      in list_init
                  for repeat    in list_repeat
                  for scaling in list_scaling
                  ]

    #systems = ['dimer_solvent_0.001']
    systems = ['for_test']

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


    for system in systems:
        system_input = InputFile( paths.get('task') + system + '/input')
        os.system(' cp -r %s/initial/ %s' % (paths.get('task') + system, paths.get('bucket')))
        initial = Dir( 'initial', paths)
        initial.checkdir()
        paths.update({'initial' :  initial.path})
        inputs.update(system_input.dict)
        for ndir in range(len(mega_list)):
            mega_list[ndir].update({ 'NDIR' : ndir,
                                     'SYSTEM' : system_input.dict.get('SYSTEM'),
                                     'PATHS_DICT' : paths,
                                     'INPUTS_DICT' : inputs
                                     })
        pool = Pool()
        pool.map( run_fssh, mega_list)
        #run_fssh(mega_list[0])







#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

def main(inputs, paths):
    print """
    """

    lol = [
        ['PROPAGATION', 'FSSH','BORN_OPPENHEIMER', 'TEST_HOP','FROZEN_HAMILTONIAN','CLASSICAL_PATH','GALILEAN'],
        ['COLLAPSE', 'T', 'F'],
        ['ANALYTICS', 'T', 'F'],
        ['FIRST_DIABAT', 1, 2],
        ['METHOD_RESCALING', 'SIMPLE', 'NACV'],
        ['METHOD_ADIAB_NACV', 'TEST', 'CONTRIBUTION','TOTAL','FAST' ],
        ['METHOD_REVERSAL', 'NEVER','ALWAYS','TRHULAR','SUBOTNIK']
    ]

    list_propagation = ['FSSH','BORN_OPPENHEIMER', 'TEST_HOP','FROZEN_HAMILTONIAN','CLASSICAL_PATH','GALILEAN']
    #list_propagation = ['FSSH']
    list_analytics   = ['T', 'F']
    #list_analytics = ['T']
    list_collapse = ['T', 'F']
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
                    'COLLAPSE'    : collapse,
                    'ANALYTICS'   : analytics,
                    'FIRST_DIABAT': diabat,
                    'METHOD_RESCALING' : rescaling,
                    'METHOD_ADIAB_NACV' : nacv,
                    'METHOD_REVERSAL'   : reversal
                   } 
                  for prop in list_propagation
                  for collapse in list_collapse 
                  for analytics in list_analytics
                  for diabat    in list_first_diabat
                  for rescaling in list_rescaling
                  for nacv      in list_nacv
                  for reversal  in list_reversal
                  ] 

    systems = ['dimer', 'trimer', 'dimer_solvent']

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



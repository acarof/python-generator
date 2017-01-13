#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

def main(inputs, paths):
    print """
    """ 

#    inputs.update({'STEPS' : 2})

    list_propagation = ['BORN_OPPENHEIMER']
    list_collapse =    ['T']
    list_first_diabat = [1]
    list_rescaling    = ['NACV']
    list_nacv         = [ 'TOTAL']
    list_reversal     = ['NEVER']
    #list_hab          = [0.0065, 0.0325, 0.0652, 0.3912, 0.9780]
    list_hab = [0.0065, 0.9780]
    list_timestep     = [0.05, 0.1, 0.5]

    mega_list = [ { 'PROPAGATION' : prop,
                    'COLLAPSE'    : collapse,
                    'METHOD_RESCALING' : rescaling,
                    'METHOD_ADIAB_NACV' : nacv,
                    'METHOD_REVERSAL'   : reversal,
                    'SCALING'           : scaling,
                    'TIMESTEP'          : timestep
                   } 
                  for prop in list_propagation
                  for collapse in list_collapse 
                  for rescaling in list_rescaling
                  for nacv      in list_nacv
                  for reversal  in list_reversal
                  for scaling   in list_hab
                  for timestep  in list_timestep
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
        #inputs.update({'STEPS': 1000})
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
                pass
            #    os.system('rm -rf run-%d' % (ndir - 1))
            else:
                print " ERROR IN CP2K FOR THOSE PARAMETERS:"
                print dict
                print " TO BE COMPARED WITH:"
                print dict_prev
                sys.exit()
            dict_prev = dict

    
    os.system('cp %sanalyser.py %s' % (paths.get('task'), paths.get('bucket')))

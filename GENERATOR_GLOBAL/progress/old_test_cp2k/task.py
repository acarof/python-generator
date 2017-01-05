#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

def main(inputs, paths):
    print """
    """ 

    inputs.update({'STEPS' : 2})

    #list_propagation = ['FSSH','BORN_OPPENHEIMER', 'TEST_HOP','FROZEN_HAMILTONIAN','CLASSICAL_PATH','GALILEAN']
    list_propagation = ['FSSH']
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

    dict_dimer = {
         'NATOMS' : 12,
         'NATOM_MOL' : 6,
         'VECTA'     : [3.527, 0.784, -0.166],
         'VECTB'     : [0, 0, 0],
         'VECTC'     : [0, 0, 0],
         'SIZE_CRYSTAL' : [2, 1, 1],
         'COORD_CHARGE' : [2, 1, 1],
         'SYSTEM': 'CRYSTAL',
         'FILE_INIT'    : 'initial_dimer',
         'MOL_NAME'     : 'ETHYLENE'
                  }
    dict_trimer = {
         'NATOMS' : 18,
         'NATOM_MOL' : 6,
         'VECTA'     : [3.527, 0.784, -0.166],
         'VECTB'     : [0, 0, 0],
         'VECTC'     : [0, 0, 0],
         'SIZE_CRYSTAL' : [3, 1, 1],
         'COORD_CHARGE' : [2, 1, 1],
         'SYSTEM': 'CRYSTAL',
         'FILE_INIT'    : 'initial_trimer',
         'MOL_NAME'     : 'ETHYLENE'
                  }
    dict_dimer_solvent = {
        'NATOMS': 135,
        'NATOM_MOL': 6,
        'SIZE_BOX'  : [30.0, 30.0, 30.0],
        'VECTA': [3.527, 0.784, -0.166],
        'VECTB': [0, 0, 0],
        'VECTC': [0, 0, 0],
        'SIZE_CRYSTAL': [2, 1, 1],
        'COORD_CHARGE': [2, 1, 1],
        'FILE_INIT': 'initial_dimer_solvent',
        'SYSTEM' : 'SOLVENT',
        'MOL_NAME': 'ETHYLENE',
        'NAME_SOLVENT' : 'AR',
        'SOLVENT'      : 'Ar'
    }

    #dict_list = [dict_dimer, dict_trimer, dict_dimer_solvent]
    dict_list = [dict_dimer_solvent]

    ndir= 0

    for structure in dict_list:
        os.system(' cp -r %s/%s %s' % (paths.get('task'), structure.get('FILE_INIT'), paths.get('bucket')))
        initial = Dir( structure.get('FILE_INIT'), paths)
        initial.checkdir()     
        paths.update({'initial' :  initial.path})
        inputs.update(structure)
        dict_prev = {}
        for dict in mega_list:
            if dict.get('PROPAGATION') == 'TEST_HOP':
               inputs.update({'STEPS':1})
            if structure.get('FILE_INIT') != 'initial_dimer':
               dict.update({'ANALYTICS' : 'F' })

            system = inputs.get('SYSTEM')
            if system == 'CRYSTAL':
                from utils import CP2KOSFSSH as Config
            elif system == 'SOLVENT':
                from utils import CP2KOSwSolventFSSH as Config
            else:
                sys.exit()

            config = Config( inputs, paths, INIT = 1, **dict)
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




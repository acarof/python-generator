#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

def main(inputs, paths):
    print """
    """ 

#    inputs.update({'STEPS' : 2})

    list_propagation = ['FSSH','BORN_OPPENHEIMER']
    list_collapse =    ['T']
    list_first_diabat = [1]
    list_rescaling    = ['NACV']
    list_nacv         = [ 'TOTAL']
    list_reversal     = ['NEVER']
    list_hab          = [0.0065, 0.0325, 0.0652, 0.3912, 0.9780]
    list_timestep     = [0.05, 0.1, 0.5, 1, 2]

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

    dict_dimer = {
         'NATOMS' : 12,
         'NATOM_MOL' : 6,
         'VECTA'     : [3.527, 0.784, -0.166],
         'VECTB'     : [0, 0, 0],
         'VECTC'     : [0, 0, 0],
         'SIZE_CRYSTAL' : [2, 1, 1],
         'COORD_CHARGE' : [2, 1, 1],
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
         'FILE_INIT'    : 'initial_trimer',
         'MOL_NAME'     : 'ETHYLENE'
                  }

    dict_list = [dict_dimer, dict_trimer]


    ndir= 0

    for structure in dict_list:
        os.system(' cp -r %s/%s %s' % ( paths.get('task'), structure.get('FILE_INIT'), paths.get('bucket')) )
        initial = Dir( structure.get('FILE_INIT'), paths)
        initial.checkdir()     
        paths.update({'initial' :  initial.path})
        inputs.update(structure)
        for dict in mega_list:
            config = CP2KFSSH( inputs, paths, INIT = 1, **dict)
            ndir = config.run(ndir)
            if ndir < 0:
               print " ERROR IN CP2K FOR THOSE PARAMETERS:"
               print dict

    
    os.system('cp %sanalyser.py %s' % (paths.get('task'), paths.get('bucket')))

#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

def main(inputs, paths):
    print """
    """

    # DEFINE PATH/EXE VARIABLES
    exe_path = '/scratch/grudorff/antoine/bin'
    paths = {'cp2k': exe_path + '/cp2k.sopt'}

    task = {
        'KIND_RUN'   : 'TASK143-ENERGETICS-FIST',
        'NUMBER_RUN' : 2,
        'NVT_STEPS'  : 7,
        'NVE_STEPS'  : 10
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

    list_timestep     = [0.1, 0.5, 1]
    mega_list = [ {
                    'TIMESTEP'          : timestep
                   }
                  for timestep  in list_timestep
                  ] 

    systems = ['dimer']

    ndir= 0
    for system in systems:
        system_input = InputFile( paths.get('task') + system + '/input')
        os.system(' cp -r %s/initial/ %s' % (paths.get('task') + system, paths.get('bucket')))
        initial = Dir( 'initial', paths)
        initial.checkdir()
        paths.update({'initial' :  initial.path})
        inputs.update(system_input.dict)

        dict_prev = {}
        for dict in mega_list:
            if dict.get('PROPAGATION') == 'TEST_HOP':
               inputs.update({'STEPS':1})
            if system_input.dict.get('FILE_INIT') != 'initial_dimer':
               dict.update({'ANALYTICS' : 'F' })

            output = Dir('output', paths)
            output.rm_mkdir()

            system = inputs.get('SYSTEM')
            if system == 'CRYSTAL':
                from utils import OSCluster as Structure
                from utils import CP2KOS as Config
            elif system == 'SOLVENT':
                from utils import OSwSolvent as Structure
                from utils import CP2KOSwSolvent as Config
            else:
                sys.exit()

            structure = Structure(inputs, paths)

            for run in range(1, inputs.get('NUMBER_RUN')+1):
                print "GO FOR RUN %d" % ndir
                config_nvt = Config(inputs, paths, ENSEMBLE='NVT', STEPS=inputs.get('NVT_STEPS')*run)
                ndir = config_nvt.run(ndir)


                config_nve = Config(inputs, paths, ENSEMBLE='NVE', STEPS=inputs.get('NVE_STEPS'), RESTART=config_nvt.ndir, **dict)
                ndir = config_nve.run(ndir)


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

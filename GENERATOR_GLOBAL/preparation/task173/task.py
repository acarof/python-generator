#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

def main(inputs, paths):
    print """
    """



    task = {
        'KIND_RUN'   : 'TASK173-DIABAT-DENSITY',
        'NCONFIG' : 1,
        'NEQ'  : 1,
        'NPROD'  : 1000,
        'TEMPLATE_FILE' : 'FIST_TEMPLATE',
        'SYSTEM': 'SOLVENT',
        'MOL_NAME'      : 'ETHYLENE',
        'FILEMOL'       : 'ethylene.coord',
        'FILECRYSTAL'   : 'ethylene_dimer.coord',
        'NAME_SOLVENT': 'NE',
        'SOLVENT': 'Ne',
        'PERIODIC' : 'XYZ',
        'SIZE_BOX'     : [40.0, 40.0, 40.0],
        'CLOSEST_DIST'  : 5,
        'VECTA'        : [3.527, 0.784, -0.166],
        'VECTB'        : [0.0, 0.0, 0.0],
        'VECTC'        : [0.0, 0.0, 0.0],
        'SIZE_CRYSTAL' : [2, 1, 1],
        'COORD_CHARGE' : [2, 1, 1],
        'RCUT'         : 12
    }
    inputs.update(task)

    list_density     = [0.03, 0.04]
    mega_list = [ {
                    'DENSITY'          : density
                   }
                  for density  in list_density
                  ]



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
    for dict in mega_list:
        inputs.update(dict)

        output = Dir('output', paths)
        output.rm_mkdir()

        print "1. CONSTRUCT THE ORGANIC CRYSTAL."
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

        print "2. RUN CP2K"
        print "GO FOR RUN %d" % ndir
        config_nvt = Config(inputs, paths, ENSEMBLE='NVT', STEPS=inputs.get('NEQ'), RESTART=None)
        ndir = config_nvt.run(ndir)

        config_nve = Config(inputs, paths, ENSEMBLE='NVE', STEPS=inputs.get('NPROD'), RESTART=config_nvt.ndir)
        ndir = config_nve.run(ndir)

        print "3. GATHER VELOCITIES AND COORDINATES, PREPARE FSSH INPUT FILES"
        fssh_parcel = FSSHParcel(inputs, paths)
        fssh_parcel.gather(config_nve.ndir)

        if os.path.exists('run-%d' % config_nve.ndir ):
            os.rename('output', 'output-%f' % inputs.get('DENSITY'))
        else:
            print " ERROR IN CP2K FOR THOSE PARAMETERS:"
            print 'run-%d' % (ndir -1)
            print dict

    os.system('cp %sanalyser.py %s' % (paths.get('task'), paths.get('bucket')))

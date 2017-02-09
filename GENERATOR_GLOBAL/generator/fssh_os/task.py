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


    os.system(' cp -r %s/%s %s' % (paths.get('task'), inputs.get('FILE_INIT'), paths.get('bucket')))
    initial = Dir(inputs.get('FILE_INIT'), paths)
    initial.checkdir()
    paths.update({'initial': initial.path})

    system = inputs.get('SYSTEM')
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    number_init = inputs.get('NUMBER_INIT', 1)
    number_random = inputs.get('NUMBER_RANDOM', 5)
    ndir= 0

    for init in range(1, number_init + 1):
        for random in range(1, number_random + 1):
            config = Config( inputs, paths, INIT = init)
            ndir = config.run(ndir)


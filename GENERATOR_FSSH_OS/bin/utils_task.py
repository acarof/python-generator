# standard modules
import string, re, struct, sys, math, os, time
import numpy as np
import imp
from collections import deque
from random import seed, randint

# custom modules
from utils import *



def generate_initial_structure(system_info, paths):
    system = system_info['SYSTEM']
    if system == 'PBC_CRYSTAL':
        from utils import OSCrystal as Structure
    else:
        raise SystemExit

    output_structure = Dir('initial/for-%s' % paths['bucket'].split('/')[-1], paths = paths, target = 'crystal')
    output_structure.rm_mkdir()

    structure = Structure(system_info, paths)
    structure.construct_organic_crystal()



def run_fist(system_info, cp2k_info, paths, steps, ndir = 0, restart_info = None, velocities = True, ensemble = 'NVE'):
    system = system_info['SYSTEM']
    if system == 'PBC_CRYSTAL':
        from utils import FISTOSCrystal as Config
    else:
        raise SystemExit

    print "GO FOR RUN %d" % ndir
    system_info.update(cp2k_info)
    config = Config(system_info, paths, ENSEMBLE=ensemble, STEPS=steps, RESTART= restart_info, VELOCITIES=velocities)
    return config.run(ndir), ndir



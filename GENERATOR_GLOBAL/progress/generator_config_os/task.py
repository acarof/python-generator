#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import Dir, FSSHParcel

def main(inputs, paths):
    print """
    START THE DIABAT OF SAMPLING OF AN ORGANIC CRYSTAL 
    """
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

    ndir = 0
    print "2. RUN CP2K"
    config_nvt = Config(inputs, paths, ENSEMBLE = 'NVT', STEPS = inputs.get('NEQ'))
    ndir = config_nvt.run(ndir)
 
    config_nve = Config(inputs, paths, ENSEMBLE = 'NVE', STEPS = inputs.get('NPROD'))
    ndir = config_nve.run(ndir)
 
 
    print "3. GATHER VELOCITIES AND COORDINATES, PREPARE FSSH INPUT FILES"
    fssh_parcel = FSSHParcel(inputs, paths)
    fssh_parcel.gather(config_nve.ndir)
 
 
    print """
    SAMPLING OF THE DIABATIC OS IS OVER!
    """

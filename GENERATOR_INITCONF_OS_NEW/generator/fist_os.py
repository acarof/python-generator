#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

try:
         finput_name = sys.argv[1]
except: 
          raise SystemExit

# DEFINE PATH/EXE VARIABLES
pexe = '/scratch/grudorff/antoine/bin'
paths = {'cp2k' : pexe + '/cp2k-test.sopt' }

# OPEN AND READ THE INPUT FILE
input_file = InputFile(finput_name)
input_file.read()

# SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
bucket = Bucket(input_file.dict)
bucket.name()
paths.update( {'bucket' : bucket.path} )

templates = Dir('templates', paths)
templates.checkdir()
templates.clean()

output = Dir('output', paths)
output.rm_mkdir()

print """
START THE DIABAT OF SAMPLING OF AN ORGANIC CRYSTAL 
""" 

print "1. CONSTRUCT THE ORGANIC CRYSTAL."
system = input_file.dict.get('SYSTEM')
if system == 'CRYSTAL':
   structure = Crystal(input_file.dict, paths)
elif system == 'SOLVENT':
   structure = Solvent(input_file.dict)
structure.create_dir()
structure.complete_dict(input_file.dict)
structure.print_info()
structure.write()

ndir = 0
print "2. RUN CP2K"
config_nvt = CP2KFistInput(input_file.dict, paths, ENSEMBLE = 'NVT', STEPS = input_file.dict.get('NEQ'))
config_nvt.print_info()
config_nvt.write()
ndir = config_nvt.run(ndir)

config_nve = CP2KFistInput(input_file.dict, paths, ENSEMBLE = 'NVE', STEPS = input_file.dict.get('NPROD'))
config_nve.print_info()
config_nve.write()
ndir = config_nve.run(ndir)


print "3. GATHER VELOCITIES AND COORDINATES, PREPARE FSSH INPUT FILES"
fssh_parcel = FSSHParcel(input_file.dict, paths)
fssh_parcel.create()
fssh_parcel.gather_vel_coord(config_nve.ndir)
fssh_parcel.gather_templates()


print """
SAMPLING OF THE DIABATIC OS IS OVER!
"""

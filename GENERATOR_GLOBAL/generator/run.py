#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

import fist_os, fssh_os, test_cp2k, energetics
from utils import *


try:
         input_name = sys.argv[1]
except:
          print "AN INPUT FILE IS REQUIRED!"
          raise SystemExit

# DEFINE PATH/EXE VARIABLES
exe_path = '/scratch/grudorff/antoine/bin'
paths = {'cp2k' : exe_path + '/cp2k-test.sopt' }

# OPEN AND READ THE INPUT FILE
input = InputFile(input_name)
input.read()

# SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
bucket = Bucket(input.dict)
bucket.name()
paths.update( {'bucket' : bucket.path} )

templates = Dir('templates', paths)
templates.checkdir()
templates.clean()

generator = Dir('generator', paths)
generator.checkdir()

kind_run = input.dict.get('KIND_RUN')
if (kind_run == 'FIST_OS'):
   output = Dir('output', paths)
   output.rm_mkdir()
   fist_os.main(input.dict, paths)   
elif (kind_run == 'FSSH_OS'):
   initial = Dir('initial', paths)
   initial.checkdir()
   fssh_os.main(input.dict, paths)   
elif (kind_run == 'TEST_CP2K'):
   test_cp2k.main(input.dict, paths)   
elif (kind_run == 'ENERGETICS'):
   energetics.main(input.dict, paths)   
else:
   print "NO METHODS TO DO!"



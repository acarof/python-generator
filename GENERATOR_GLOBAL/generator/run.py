#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import importlib, imp

import fist_os, fssh_os, test_cp2k
from utils import *


try:
         input_info = sys.argv[1]
except:
          print "AN INPUT FILE IS REQUIRED!"
          raise SystemExit

# DEFINE PATH/EXE VARIABLES
exe_path = '/scratch/grudorff/antoine/bin'
paths = {'cp2k' : exe_path + '/cp2k-test.sopt' }

# OPEN AND READ THE INPUT FILE
if 'task/' in input_info:
   input = InputFile(input_info + '/input')
   input.read()
   kind_run = 'TASK'
else:
   input = InputFile(input_name)
   input.read()
   kind_run = input.dict.get('KIND_RUN', 'NO_METHOD')

# SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
bucket = Bucket(input.dict)
bucket.name()
paths.update( {'bucket' : bucket.path} )

templates = Dir('templates', paths)
templates.checkdir()
templates.clean()

generator = Dir('generator', paths)
generator.checkdir()

do_clean = input.dict.get('CLEAN')
if do_clean == 'YES':
   os.system('rm -rf run-*' )
   os.system('rm -rf fail-*' )


if (kind_run == 'NO_METHOD'):
   print "A METHOD IS REQUIRED!"
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
elif (kind_run == 'TASK'):
   task = Dir(input_info)
   paths.update( {'task' : task.path} )
   #mod = imp.load_source(   hs.get('generator') + kind_run.lower() + '.py')
   print task.path + 'task.py'
   mod = imp.load_source('task', task.path + 'task.py')
   mod.main(input.dict, paths)
else:
   print "METHOD %s IS UNKNOWN" % kind_run


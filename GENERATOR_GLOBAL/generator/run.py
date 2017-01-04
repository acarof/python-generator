#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import importlib, imp


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
if input_info == 'clean':
   os.system('rm -rf run-*')
   os.system('rm -rf fail-*')
   os.system('rm -rf output')
   os.system('rm -rf initial_*')
   os.system('rm -f analyser.py')
   os.system('rm -rf tmp')
   sys.exit()
elif 'task/' in input_info:
   input = InputFile(input_info + '/input')
   kind_run = 'TASK'
elif 'progress/' in input_info:
   input = InputFile(input_info + '/input')
   input.dict.update({'TEST' : 'YES'})
   kind_run = 'TASK'
else:
   input = InputFile(input_name)
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
   mod = imp.load_source('task', task.path + 'task.py')
   mod.main(input.dict, paths)
else:
   print "METHOD %s IS UNKNOWN" % kind_run


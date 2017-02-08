#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import importlib, imp


from utils import *


try:
         input_info = sys.argv[1]
except:
          print "A INPUT DIRECTORY IS REQUIRED!"
          raise SystemExit

# OPEN AND READ THE INPUT FILE
if 'task/' in input_info:
   if os.path.isfile(input_info + '/input'):
      input = InputFile(input_info + '/input', True)
   else:
      input = InputFile(input_info + '/input', False)
   input.dict.update({'TEST': 'NO'})
   kind_run = 'TASK'
elif 'progress/' in input_info:
   if os.path.isfile(input_info + '/input'):
      input = InputFile(input_info + '/input', True)
   else:
      input = InputFile(input_info + '/input', False)
   input.dict.update({'TEST' : 'YES'})
   kind_run = 'TASK'
else:
   input = InputFile(input_name)
   kind_run = input.dict.get('KIND_RUN', 'NO_METHOD')

# SELECT THE METHOD TO RUN
paths = {}
if (kind_run == 'NO_METHOD'):
   print "A METHOD IS REQUIRED!"
elif (kind_run == 'FIST_OS'):
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
   input.dict.update({'INPUT_INFO': input_info})
   mod = imp.load_source('task', task.path + 'task.py')
   mod.main(input.dict, paths)
else:
   print "METHOD %s IS UNKNOWN" % kind_run


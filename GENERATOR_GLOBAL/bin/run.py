#!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy
import importlib, imp

# custom modules
from utils import *
import machine



try:
    input_info = sys.argv[1]
except:
    print "A INPUT DIRECTORY IS REQUIRED!"
    raise SystemExit

try:
    nworker = int( sys.argv[2] )
except:
    nworker = -1

# PATHS CONTAINS ALL THE PATHS
paths = {}


# OPEN AND READ THE INPUT FILE
def read_input(path):
   if os.path.isfile(path + '/input'):
      input = InputFile(path + '/input')
   else:
      input = InputFile('NONE')
   return input


if 'generator/' in input_info:
   input = read_input(input_info)
   input.dict.update({'TEST' : 'YES'})
elif 'preparation/' in input_info:
   input = read_input(input_info)
   input.dict.update({'TEST' : 'YES'})
else:
   input = read_input(input_info)
   input.dict.update({'TEST': 'NO'})
input.dict.update({'INPUT_INFO': input_info})
input.dict.update({'NWORKER' : nworker})

# GET CP2K READY
paths.update(
   { 'cp2k' : machine.get_cp2k_path()
     })
#machine.source_cp2k()


# UPLOAD task.py as a MODULE
task = Dir(input_info)
paths.update( {'task' : task.path} )
mod = imp.load_source('task', task.path + 'task.py')


# RUN MAIN FUNCTION OF task.py
mod.main(input.dict, paths)


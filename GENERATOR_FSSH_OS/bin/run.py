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
    print "A TASK.py IS REQUIRED!"
    raise SystemExit

try:
    nworker = int( sys.argv[2] )
except:
    nworker = -1

try:
    ncore = int( sys.argv[3])
except:
    ncore = -1

# OPEN AND READ THE INPUT FILE
input = InputFile('NONE')
input.dict.update({'INPUT_INFO': input_info})
input.dict.update({'NWORKER' : nworker})
input.dict.update({'NCORE' : ncore})

paths = {}
# GET CP2K READY, get the correct machine to use by bin/machine.py
paths.update(
   { 'cp2k' : machine.get_cp2k_path()
     })
#machine.source_cp2k()


# UPLOAD task.py as a MODULE
task = Dir(input_info)
paths.update( {'task' : task.path} )
mod = imp.load_source('tasks', input_info)


# RUN MAIN FUNCTION OF task.py
mod.main(input.dict, paths)


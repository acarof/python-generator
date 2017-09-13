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
def create_cp2k_path():
    print "The file cp2k.path doesn't exist, please provide the path for CP2K:"
    path = raw_input('> ')
    cp2k_path = open('bin/cp2k.path', 'w')
    cp2k_path.write(path)
    cp2k_path.close()
    return open('bin/cp2k.path', 'r')
try:
    cp2k_path = open('bin/cp2k.path', 'r')
except:
    cp2k_path = create_cp2k_path()
paths.update(
    {'cp2k' : cp2k_path.readline()}
)


# UPLOAD task.py as a MODULE
task = Dir(input_info)
paths.update( {'task' : task.path} )
mod = imp.load_source('tasks', input_info)


# RUN MAIN FUNCTION OF task.py
mod.main(input.dict, paths)


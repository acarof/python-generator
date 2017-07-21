#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import importlib, imp


from utils import *

name_bucket = os.getcwd().split('/')[-1]
if "GENERATOR" not in os.getcwd():
    list = name_bucket.split('-')
    list[-1] = 'CLEAN'
    new_name =  '-'.join(list)
    os.chdir('..')
    os.system('mv %s %s' % (name_bucket, new_name))
    os.chdir(new_name)

os.system('rm -rf run-*')
os.system('rm -rf fail-*')
os.system('rm -rf output*')
#os.system('rm -rf initial*')
os.system('rm -f analyser.py')
os.system('rm -rf tmp*')
os.system('rm -rf iagodb.json')
os.system('rm -rf None')
os.system('rm -rf info.*')



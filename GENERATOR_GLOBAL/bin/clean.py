#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy
import importlib, imp


from utils import *


os.system('rm -rf run-*')
os.system('rm -rf fail-*')
os.system('rm -rf output')
os.system('rm -rf initial_*')
os.system('rm -f analyser.py')
os.system('rm -rf tmp')
os.system('rm -rf iagodb.json')



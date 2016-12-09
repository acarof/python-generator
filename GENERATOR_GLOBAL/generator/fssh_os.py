#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

def main(inputs, paths):
    print """
    """ 

    number_init = inputs.get('NUMBER_INIT', 1)
    number_random = inputs.get('NUMBER_RANDOM', 5)
    ndir= 0

    for init in range(1, number_init + 1):
        for random in range(1, number_random + 1):
            config = CP2KFSSH( inputs, paths, INIT = init)
            ndir = config.run(ndir)


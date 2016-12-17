#!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

from utils import *

def main(inputs, paths):
    print """
    """

    os.system(' cp -r %s/%s %s' % (paths.get('task'), inputs.get('FILE_INIT'), paths.get('bucket')))
    initial = Dir(inputs.get('FILE_INIT'), paths)
    initial.checkdir()
    paths.update({'initial': initial.path})


    number_init = inputs.get('NUMBER_INIT', 1)
    number_random = inputs.get('NUMBER_RANDOM', 5)
    ndir= 0

    for init in range(1, number_init + 1):
        for random in range(1, number_random + 1):
            config = CP2KOSFSSH( inputs, paths, INIT = init)
            ndir = config.run(ndir)


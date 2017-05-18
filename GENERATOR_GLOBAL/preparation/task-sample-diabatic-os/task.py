#!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy

# custom modules
from utils import Dir, FSSHParcel, Bucket
from utils_task import run_fist_nvt_nve_extract




def main(inputs, paths):
    print """
    START THE DIABAT OF SAMPLING OF AN ORGANIC CRYSTAL 
    """

    # SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
    bucket = Bucket(inputs)
    bucket.name()
    paths.update({'bucket': bucket.path})

    templates = Dir('templates', paths)
    templates.checkdir()
    templates.clean()

    bin = Dir('bin', paths)
    bin.checkdir()

    output = Dir('output', paths)
    output.rm_mkdir()


    # PREPARE THE MEGA_LIST FOR POOL
    mega_list = [{}]
    for ndir in range(len(mega_list)):
        mega_list[ndir].update({ 'NDIR' : ndir*2,
                                 'PATHS_DICT' : paths,
                                 'INPUTS_DICT' : inputs
                                 })


    # RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
    nworker = 0
    #nworker = inputs['NWORKER']
    #if nworker == -1:
    #    nworker = cpu_count()

    if nworker == 0:
        for dict_ in mega_list:
            run_fist_nvt_nve_extract(dict_)
    #else:
    #    pool = Pool(nworker)
    #     pool.map(run_fist, mega_list)


    # COPY ANALYSER
    os.system('cp %sanalyser.py %s' % (paths.get('task'), paths.get('bucket')))
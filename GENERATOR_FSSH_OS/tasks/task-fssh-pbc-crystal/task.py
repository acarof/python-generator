#!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy as np
from multiprocessing import Pool, cpu_count
import imp
from itertools import product as iterprod


# custom modules
from utils import *
from utils_task import *
from find_crystal_bb import find_crystal_bb


def main(task_info, paths):
    print """
    """

    task = {
        'KIND_RUN' : 'TONAME',
        'FILE_INIT': 'GENERATOR_FSSH_OS',
        'LENGTH_FS': 1,
        'INITIALIZATION': 'DIABATIC',
        'NUMBER_CONFIG'        : 1,
        'NUMBER_REPEAT'  :  1,
        'LIGHT' : False,
        #'LIST_ACTIVATED' : [1,2,3]
    }
    task_info.update(task)

    system_info = (InputFile('%s/system.info' % 'initial/from-%s' % task_info['FILE_INIT']).dict)
    system_info.update({
        'AOM_RADIUS' : 3.0
    })

    cp2k_param = [
        #[ 'PROPAGATION', 'FSSH', 'BORN_OPPENHEIMER'],
        ['PROPAGATION', 'FSSH'],
        #['DECO', 'NO_DECO_CORR','INSTANT_COLLAPSE','DAMPING']
        [ 'DECO', 'DAMPING'],
        #['METHOD_RESCALING','NACV','SIMPLE_QSYS']
        [ 'METHOD_RESCALING', 'NACV'],
        #[ 'METHOD_ADIAB_NACV', 'FAST','TOTAL']
        [ 'METHOD_ADIAB_NACV', 'FAST'],
        #['METHOD_REVERSAL', 'NEVER', 'ALWAYS', 'TRUHLAR', 'SUBOTNIK']
        [ 'METHOD_REVERSAL', 'ALWAYS'],
        #[ 'SCALING', 0.05, 0.03, 0.01, 0.008, 0.005, 0.003, 0.001, 0.0005, 0.0001, 0.00005],
        [ 'SCALING', 0.00005],
        #['SCALING', 0.003, 0.00005],
        #[ 'TIMESTEP', 0.01, 0.05, 0.1, 0.5],
        [ 'TIMESTEP', 0.1],
        #['EDC_E0', 0.01, 0.1, 1.0],
        ['EDC_E0', 0.1],
        #['ELECTRONIC_STEPS', 5, 10, 50],
        ['ELECTRONIC_STEPS', 5],
        [ 'TEMPLATE_FILE', 'FSSH_PBC_CRYSTAL.template'],
        ['FORCEFIELD_FILE', 'ANTRACENE_FF.prm'],
        ['INITIALIZATION', 'DIABATIC'],
        ['INIT', 0],
        ['REPEAT'] + range(task_info.get('NUMBER_REPEAT'))
    ]


    # FIND ACTIVE MOLECULES WITH GUIDO'S TOOL
    system_info.update({
        'LIST_ACTIVATED' : find_crystal_bb(
            abc=system_info['ABC'] + system_info['ALPHA_BETA_GAMMA'],
            start= system_info['STARTING_POINT'],
            length = system_info['LENGTH'],
            vector = system_info['DIRECTION'],
            radius = system_info['RCUT'],
            has_atoms = True,
            radius_aom = system_info['AOM_RADIUS'],
            psf_file = '%s/input-1.psf' % 'initial/from-%s' % task_info['FILE_INIT'],
            xyz_file = '%s/crystal.xyz' % 'initial/from-%s' % task_info['FILE_INIT']

        )
    })

    print system_info['LIST_ACTIVATED']
    #system_info['LIST_ACTIVATED'] = [1, 2, 3]
    #exit()

    # BUILD THE MEGA_LISTS
    second_list = [ sublist[1:] for sublist in cp2k_param]
    total_list  = list(iterprod(*second_list))
    mega_list = []
    for sublist in total_list:
        subdict = {}
        for index in range(len(sublist)):
            subdict.update({
                cp2k_param[index][0] : sublist[index]
            })
        mega_list.append(subdict)


    # SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
    bucket = Bucket(task_info)
    bucket.name()
    paths.update({'bucket': bucket.path})

    task = Dir(task_info.get('INPUT_INFO'))
    paths.update( {'task' : task.path} )


    # PREPARE THE MEGA_LIST FOR POOL
    for ndir in range(len(mega_list)):
        mega_list[ndir].update({ 'NDIR' : ndir,
                                 'PATHS_DICT' : paths,
                                 'INPUTS_DICT' : task_info
                                 })


    # RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
    nworker = task_info['NWORKER']
    if nworker == 0:
        for cp2k_info in mega_list:
            cp2k_info.update(system_info)
            run_fssh(cp2k_info)
    else:
        from multiprocessing import Pool, cpu_count
        if nworker == -1:
            nworker = cpu_count()
        pool = Pool(nworker)
        pool.map( run_fssh, mega_list)









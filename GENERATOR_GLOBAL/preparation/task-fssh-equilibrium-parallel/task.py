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



def main(task_info, paths):
    print """
    """

    task = {
        'KIND_RUN' : 'TONAME',
        'FILE_INIT': 'TASK271-SAMPLE-BO-30ps-170517-b7905f4556ed79fc2323eeec5105416e',
        'LENGTH_FS': 10,
        'INITIALIZATION': 'ADIABATIC',
        'NUMBER_ADIABAT' : 2,
        'NUMBER_CONFIG'        : 1,
        'NUMBER_REPEAT'  :  5,
        'LIGHT' : True
    }
    task_info.update(task)


    cp2k_param = [
        [ 'PROPAGATION', 'FSSH'],
        #['DECO', 'NO_DECO_CORR','INSTANT_COLLAPSE','DAMPING']
        [ 'DECO', 'DAMPING'],
        #['METHOD_RESCALING','NACV','SIMPLE_QSYS']
        [ 'METHOD_RESCALING', 'NACV'],
        #[ 'METHOD_ADIAB_NACV', 'FAST','TOTAL']
        [ 'METHOD_ADIAB_NACV', 'FAST'],
        #['METHOD_REVERSAL', 'NEVER', 'ALWAYS', 'TRUHLAR', 'SUBOTNIK']
        [ 'METHOD_REVERSAL', 'ALWAYS'],
        #[ 'SCALING', 0.05, 0.03, 0.01, 0.008, 0.005, 0.003, 0.001, 0.0005, 0.0001, 0.00005],
        [ 'SCALING', 0.03, 0.02, 0.01, 0.008, 0.003, 0.0005, 0.00005],
        #[ 'TIMESTEP', 0.01, 0.05, 0.1, 0.5],
        [ 'TIMESTEP', 0.1],
        #['EDC_E0', 0.01, 0.1, 1.0],
        ['EDC_E0', 0.1],
        #['ELECTRONIC_STEPS', 5, 10, 50],
        ['ELECTRONIC_STEPS', 5],
        [ 'TEMPLATE_FILE', 'FSSH_CORE.template'],
        [ 'FORCEFIELD_FILE', 'FSSH_FF.template'],
        ['INITIALIZATION', 'ADIABATIC'],
        ['INIT'] + range(1, task_info.get('NUMBER_CONFIG') + 1),
        ['REPEAT'] + range(task_info.get('NUMBER_REPEAT'))
    ]

    # This data are taken from: extract-scaling-adiabat-TASK271-SAMPLE-BO-30ps-170517-b7905f4556ed79fc2323eeec5105416e-1705172113
    dict_for_equilibrium = { \
        5e-05: [0.882946564262, 0.117053435738], \
        0.0001: [0.883844010925, 0.116155989075], \
        0.0005: [0.89086291711, 0.10913708289], \
        0.001: [0.899233079636, 0.100766920364], \
        0.003: [0.928251700237, 0.0717482997627], \
        0.005: [0.950477468394, 0.0495225316059], \
        0.008: [0.973033713426, 0.0269662865739], \
        0.01: [0.98253952535, 0.0174604746505], \
        0.02: [0.998429468654, 0.00157053134578], \
        0.03: [0.999885479206, 0.000114520794399], \
        0.05: [0.999999514713, 4.85286930795e-07]
    }


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

    templates = Dir('templates', paths)
    templates.checkdir()
    templates.clean()

    supinitial = Dir('initial')
    supinitial.checkdir()

    bin = Dir('bin', paths)
    bin.checkdir()


    # PREPARE THE MEGA_LIST FOR POOL
    for ndir in range(len(mega_list)):
        mega_list[ndir].update({ 'NDIR' : ndir,
                                 'PATHS_DICT' : paths,
                                 'INPUTS_DICT' : task_info,
                                 'DICT_EQUILIBRIUM' : dict_for_equilibrium
                                 })


    # RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
    nworker = task_info['NWORKER']
    if nworker == 0:
        for cp2k_info in mega_list:
            run_fssh(cp2k_info)
    else:
        from multiprocessing import Pool, cpu_count
        if nworker == -1:
            nworker = cpu_count()
        pool = Pool(nworker)
        pool.map( run_fssh, mega_list)









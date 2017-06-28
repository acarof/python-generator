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
        'FILE_INIT': 'TASK279-SAMPLE-BO-CORRECT-TEMP-100ps-20dabab-170601-8da556911f813378f0577abfd206e148',
        'FILE_DICT' : 'TASK279-SAMPLE-BO-CORRECT-TEMP-100ps-20dabab-170601-8da556911f813378f0577abfd206e148-1706121812',
        'LENGTH_FS': 10000,
        'INITIALIZATION': 'ADIABATIC',
        'NUMBER_ADIABAT' : 2,
        'NUMBER_CONFIG'        : 1,
        'FIRST_CONFIG' : 200,
        'FINAL_CONFIG' : 500,
        'NUMBER_REPEAT'  :  1,
        'LIGHT' : True
    }
    task_info.update(task)

    seed()
    cp2k_param = [
        [ 'PROPAGATION', 'FSSH'],
        #['SURFACE_HOP_CHOICE', 'TRIVIAL_HOP_CORRECT','UNMODIFIED_SURF_HOP']
        ['SURFACE_HOP_CHOICE', 'TRIVIAL_HOP_CORRECT'],
        # ['DECO', 'NO_DECO_CORR','INSTANT_COLLAPSE','DAMPING']
        [ 'DECO', 'DAMPING'],
        #['METHOD_RESCALING','NACV','SIMPLE_QSYS']
        [ 'METHOD_RESCALING', 'NACV'],
        #[ 'METHOD_ADIAB_NACV', 'FAST','TOTAL']
        [ 'METHOD_ADIAB_NACV', 'FAST'],
        #['METHOD_REVERSAL', 'NEVER', 'ALWAYS', 'TRUHLAR', 'SUBOTNIK']
        [ 'METHOD_REVERSAL', 'ALWAYS'],
        #[ 'SCALING', 0.05, 0.03, 0.01, 0.008, 0.005, 0.003, 0.001, 0.0005, 0.0001, 0.00005],
        #[ 'SCALING', 0.03, 0.02, 0.01, 0.008, 0.003, 0.0005, 0.00005],
        ['SCALING', 0.03],
        #['SCALING', 0.003, 0.00005],
        #[ 'TIMESTEP', 0.01, 0.05, 0.1, 0.5],
        [ 'TIMESTEP', 0.5],
        #['EDC_E0', 0.01, 0.1, 1.0],
        ['EDC_E0', 0.1],
        #['ELECTRONIC_STEPS', 5, 10, 50],
        ['ELECTRONIC_STEPS', 5],
        [ 'TEMPLATE_FILE', 'FSSH_CORE.template'],
        [ 'FORCEFIELD_FILE', 'FSSH_FF.template'],
        ['INITIALIZATION', 'ADIABATIC'],
        ['INIT_CONFIG'] + [ (ind + 1, range(task_info['FIRST_CONFIG'], task_info['FINAL_CONFIG'], (task_info['FINAL_CONFIG'] - task_info['FIRST_CONFIG']) / task_info['NUMBER_CONFIG'])[ind])
                      for ind in range(task_info['NUMBER_CONFIG'])],
        ['REPEAT'] + range(task_info.get('NUMBER_REPEAT')),
    ]

    print "How many run?", len(cp2k_param)

    dict_for_equilibrium = {}
    with open('initial/result-thermodynamics-population-extract-scaling-adiabat-%s/Simulation-1.dat'  %
              task_info['FILE_DICT']) as file_dict:
        for line in file_dict.readlines()[1:]:
            if '#' in line:
                pass
            else:
                dict_for_equilibrium[
                    float(line.split()[0])
                ] = [ float(x) for x in line.split()[1:] ]

    print dict_for_equilibrium

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









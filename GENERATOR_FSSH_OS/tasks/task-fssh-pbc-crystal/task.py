#!/usr/bin/python

# standard modules
#import string, re, struct, sys, math, os, time
import numpy
#import imp
from itertools import product as iterprod


# custom modules
from utils import *
from utils_task import *
from find_crystal_bb import find_molecules


def main(task_info, paths):
    print """
    """

    task = {
        #################### CAN BE CHANGED ###############################################
        'KIND_RUN' : 'TONAME',                      # NAME OF YOUR RUN
        'FILE_INIT': 'GENERATOR_FSSH_OS-temp-100',           # NAME OF THE RUN OF INITIALIZATION
        'NCONFIG_INIT' : 10,
        'NPROD_INIT'   : 10,
        'NUMBER_CONFIG': 10,
        'NUMBER_REPEAT': 5,
        'LENGTH_FS': 20,                             # LENGTH IN FS
        ###################################################################################
        'INITIALIZATION': 'DIABATIC',
        'LIGHT': True,
    }
    task_info.update(task)


    list_config_init = range(0, task['NPROD_INIT'], task['NPROD_INIT'] / task['NCONFIG_INIT'])

    cp2k_param = [
        #################### CAN BE CHANGED ###############################################
        ['PROPAGATION', 'FSSH'],    # METHOD OF PROPAGATION: FSSH OR BORN_OPPENEHIMER
        ['SCALING', 0.06685],       # SCALING FACTOR IN HARTREE (C = 1.819 eV)
        ['TIMESTEP', 0.5],          # TIMESTEP IN FS
        ['REPEAT'] + range(task_info.get('NUMBER_REPEAT')), # NUMBER OF FSSH RUN PER STARTING POINT
         ###################################################################################
        ['INIT'] + \
        [ list_config_init[x] for x in
         range(0, len(list_config_init), len(list_config_init) / task['NUMBER_CONFIG'])], # STARTING POINTS
        [ 'DECO', 'INSTANT_COLLAPSE'],
        [ 'METHOD_RESCALING', 'NACV'],
        [ 'METHOD_ADIAB_NACV', 'FAST'],
        [ 'METHOD_REVERSAL', 'ALWAYS'],
        ['EDC_E0', 0.1],
        ['ELECTRONIC_STEPS', 5],
        [ 'TEMPLATE_FILE', 'FSSH_PBC_CRYSTAL.template'],
        ['FORCEFIELD_FILE', 'ANTRACENE_FF.prm'],
        ['INITIALIZATION', 'DIABATIC'],
    ]


    system_info = (InputFile('%s/system.info' % 'initial/from-%s' % task_info['FILE_INIT']).dict)
    system_info.update({
        'AOM_RADIUS' : 3.0
    })


    # FIND ACTIVE MOLECULES WITH GUIDO'S TOOL
    length = 0.0
    for x,y in zip(system_info['ABC'], system_info['DIRECTION']):
        length += (x*y)**2
    length= numpy.sqrt(length)*(system_info['NUMBER_MOL_ACTIVE'] - 1)
    system_info.update({
        'LIST_ACTIVATED' : find_molecules(
            coord_charge = system_info['COORD_CHARGE'],
            size_crystal = system_info['SIZE_CRYSTAL'],
            length = length,
            vector = system_info['DIRECTION'],
            radius_aom = system_info['AOM_RADIUS'],
            psf_file = '%s/input-1.psf' % 'initial/from-%s' % task_info['FILE_INIT'],
            xyz_file = '%s/crystal.xyz' % 'initial/from-%s' % task_info['FILE_INIT'],
            nmol_unit= system_info['NMOL_UNIT']

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
    paths.update({'bucket': os.getcwd()})
    task = Dir(task_info.get('INPUT_INFO'))
    paths.update( {'task' : task.path} )
    # PATHS CONTAINS ALL THE PATHS
    for directory in ['bin', 'initial', 'scripts', 'structures', 'tasks', 'templates', 'tools', 'topologies']:
        dir = Dir(directory, paths)
        dir.checkdir()

    # PREPARE THE MEGA_LIST FOR POOL
    for ndir in range(len(mega_list)):
        mega_list[ndir].update({ 'NDIR' : ndir,
                                 'PATHS_DICT' : paths,
                                 'INPUTS_DICT' : task_info
                                 })
        mega_list[ndir].update(system_info)


    # RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
    nworker = task_info['NWORKER']
    if nworker == 0:
        for cp2k_info in mega_list:
            #cp2k_info.update(system_info)
            run_fssh(cp2k_info)
    else:
        from multiprocessing import Pool, cpu_count
        if nworker == -1:
            nworker = cpu_count()
        pool = Pool(nworker)
        pool.map( run_fssh, mega_list)









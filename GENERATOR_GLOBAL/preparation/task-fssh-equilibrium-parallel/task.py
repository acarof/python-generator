#!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy as np
from multiprocessing import Pool, cpu_count
import imp
from itertools import product as iterprod


# custom modules
from utils import *



def run_fssh( cp2k_info):
    task_info = cp2k_info.get('INPUTS_DICT')
    paths  = cp2k_info.get('PATHS_DICT')


    # DETERMINE INITIAL STATE TO GET BOLTZMANN RATIO
    scaling = cp2k_info['SCALING']
    pop_ground_state = estimate_boltzmann_ratio(scaling, paths, task_info.get('NUMBER_CONFIG')*task_info['NUMBER_REPEAT'])
    indice = (cp2k_info['REPEAT'] + 1)  + (cp2k_info['INIT'] -1 )*task_info['NUMBER_REPEAT']
    if indice < pop_ground_state:
        initial = Dir('initial/state-1-scaling-%s-%s' % (scaling, task_info['FILE_INIT']), paths)
        cp2k_info.update({'FIRST_ADIABAT' : 1})
    else:
        initial = Dir('initial/state-2-scaling-%s-%s' % (scaling, task_info['FILE_INIT']), paths)
        cp2k_info.update({'FIRST_ADIABAT' : 2})


    cp2k_info.update({'STEPS' : int(task_info['LENGTH_FS'] / cp2k_info['TIMESTEP'] ) } )
    cp2k_info.update({'PRINT':  int(task_info['LENGTH_FS'] / cp2k_info['TIMESTEP']) } )
    cp2k_info.update({'PRINT_FSSH': int( 1 / cp2k_info['TIMESTEP']) })


    initial.checkdir()
    paths.update({'initial': initial.path})
    systems = InputFile( initial.path + '/system.info').dict
    cp2k_info.update(systems)


    system = cp2k_info['SYSTEM']
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()


    config = Config(cp2k_info, paths)
    ndir = cp2k_info['NDIR']
    print "GO FOR RUN %d" % ndir
    config.run(ndir)


def estimate_boltzmann_ratio(scaling, paths, nconfig):  # scaling in Ha
    jacob_coupling = {
        89: 8,
        133: 11,
        177: 15,
        266: 23,
        355: 30,
        532: 45,
        710: 61,
        887: 76,
        1330: 114,
        1774: 151,
        2217: 189,
        2661: 227,
        3104: 265,
        3548: 303,
        3991: 341,
        4435: 378,
        4878: 418,
        5322: 455
    }

    scalings = np.array(jacob_coupling.keys()) * 0.0000367493  # convert in Ha
    couplings = np.array(jacob_coupling.values()) * 0.001  # convert in eV
    fit = np.polyfit(scalings, couplings, 1)
    coupling =  fit[0] * scaling + fit[1]

    marcus = imp.load_source('scripts', paths.get('bucket') + '/scripts/marcus.py')
    temperature = 300
    reorga = 0.100
    free_energy = 0

    return np.rint( nconfig * marcus.calculate_boltzman_ratio(reorga, free_energy, coupling, temperature, state = 'Ground') )

def estimate_scaling(coupling):  # coupling in eV
    jacob_coupling = {
        89: 8,
        133: 11,
        177: 15,
        266: 23,
        355: 30,
        532: 45,
        710: 61,
        887: 76,
        1330: 114,
        1774: 151,
        2217: 189,
        2661: 227,
        3104: 265,
        3548: 303,
        3991: 341,
        4435: 378,
        4878: 418,
        5322: 455
    }

    scalings = np.array(jacob_coupling.keys()) * 0.0000367493  # convert in Ha
    couplings = np.array(jacob_coupling.values()) * 0.001  # convert in eV
    fit = np.polyfit(scalings, couplings, 1)
    scaling = (coupling - fit[1]) / fit[0]
    #coupling =  fit[0] * scaling + fit[1]

    return scaling

def round_to_1(x):
    return round(x, -int(np.floor(np.log10(np.abs(x)))))

def get_list_scaling(number, reorga):
    division = { 1 : 1, 2 : 2, 3 : 3, 4 :5, 5 : 10, 6 : 100, 7 : 1000}
    scalings = []
    for i in range(1, number + 1):
        div = division[i]
        coupling = reorga / div
        scaling = estimate_scaling(coupling)
        scalings.append(round_to_1(scaling))
    return scalings


def main(task_info, paths):
    print """
    """

    task = {
        'KIND_RUN' : 'TONAME',
        'FILE_INIT': 'TASK234-SAMPLE-TWO-ADIABATS-170405-1655c3bacaee42ccabccff93dbff0f92',
        'TEMPLATE_FILE': 'FSSH_CORE.template',
        'FORCEFIELD_FILE': 'FSSH_FF.template',
        'LENGTH_FS': 1,
        'INITIALIZATION': 'ADIABATIC',
        'NUMBER_CONFIG'        : 1,
        'NUMBER_REPEAT'  :  10,
        'NUMBER_SCALING' : 1,
        'REORGANIZATION_ENERGY' : 0.1  # eV
    }
    task_info.update(task)


    cp2k_param = [
        [ 'PROPAGATION', 'FSSH'],
        [ 'DECO', 'DAMPING'],
        [ 'METHOD_RESCALING', 'NACV'],
        [ 'METHOD_ADIAB_NACV', 'FAST'],
        [ 'METHOD_REVERSAL', 'ALWAYS'],
        [ 'INIT'] + range(1, task_info.get('NUMBER_CONFIG') + 1),
        [ 'REPEAT'] + range(task_info.get('NUMBER_REPEAT')),
        [ 'SCALING'] + get_list_scaling( task_info['NUMBER_SCALING'], task_info['REORGANIZATION_ENERGY']  ),
        [ 'TIMESTEP', 0.5],
        [ 'TEMPLATE_FILE', 'FSSH_CORE.template'],
        [ 'FORCEFIELD_FILE', 'FSSH_FF.template']
    ]


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
                                 'INPUTS_DICT' : task_info
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









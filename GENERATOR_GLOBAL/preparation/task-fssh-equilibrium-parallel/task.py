#!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy as np
from multiprocessing import Pool, cpu_count
import imp

# custom modules
from utils import *



def run_fssh( dict_):
    inputs = dict_.get('INPUTS_DICT')
    paths  = dict_.get('PATHS_DICT')

    # DETERMINE INITIAL STATE TO GET BOLTZMANN RATIO
    scaling = dict_['SCALING']
    pop_ground_state = estimate_boltzmann_ratio(scaling, paths, inputs.get('NUMBER_CONFIG')*inputs['NUMBER_REPEAT'])
    indice = (dict_['REPEAT'] + 1)  + (dict_['INIT'] -1 )*inputs['NUMBER_REPEAT']
    if indice < pop_ground_state:
        initial = Dir('initial/state-1-scaling-%s-%s' % (scaling, inputs['FILE_INIT']), paths)
        inputs.update({'FIRST_ADIABAT' : 1})
    else:
        initial = Dir('initial/state-2-scaling-%s-%s' % (scaling, inputs['FILE_INIT']), paths)
        inputs.update({'FIRST_ADIABAT' : 2})


    initial.checkdir()
    paths.update({'initial': initial.path})
    systems = InputFile( initial.path + '/system.info').dict
    inputs.update(systems)

    system = inputs['SYSTEM']
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    config = Config(inputs, paths, **dict_)
    ndir = dict_['NDIR']
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


def main(inputs, paths):
    print """
    """

    task = {
        'KIND_RUN' : 'TONAME',
        'FILE_INIT': 'TASK234-SAMPLE-TWO-ADIABATS-170405-1655c3bacaee42ccabccff93dbff0f92',
        'TEMPLATE_FILE': 'FSSH_CORE.template',
        'FORCEFIELD_FILE': 'FSSH_FF.template',
        'STEPS'    : 10000,
        'PRINT'    : 10000,
        'TIMESTEP' : 0.5,
        'PRINT_FSSH' : 1,
        'INITIALIZATION': 'ADIABATIC',
        'NUMBER_CONFIG'        : 1,
        'NUMBER_REPEAT'  :  1,
        'NUMBER_SCALING' : 1,
        'REORGANIZATION_ENERGY' : 0.1  # eV
    }
    inputs.update(task)

    list_propagation = ['FSSH']
    #list_decoherences = ['INSTANT_COLLAPSE', 'DAMPING']
    list_decoherences = ['DAMPING']
    #list_rescaling    = ['NACV','SIMPLE_QSYS']
    list_rescaling = ['NACV']
    list_nacv         = ['FAST']
    #list_reversal = ['NEVER', 'ALWAYS', 'TRHULAR', 'SUBOTNIK']
    list_reversal = ['NEVER']
    list_init     = range(1, inputs.get('NUMBER_CONFIG') + 1)
    list_repeat   = range(inputs.get('NUMBER_REPEAT'))
    list_scaling = get_list_scaling( inputs['NUMBER_SCALING'], inputs['REORGANIZATION_ENERGY']  )
    #list_scaling = [0.03]


    mega_list = [ { 'PROPAGATION' : prop,
                    'DECO': deco,
                    'METHOD_RESCALING' : rescaling,
                    'METHOD_ADIAB_NACV' : nacv,
                    'METHOD_REVERSAL'   : reversal,
                    'INIT'              : init ,
                    'REPEAT'             : repeat,
                    'SCALING'            : scaling
                   }
                  for prop in list_propagation
                  for deco in list_decoherences
                  for rescaling in list_rescaling
                  for nacv      in list_nacv
                  for reversal  in list_reversal
                  for init      in list_init
                  for repeat    in list_repeat
                  for scaling in list_scaling
                  ]



    # SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
    bucket = Bucket(inputs)
    bucket.name()
    paths.update({'bucket': bucket.path})

    task = Dir(inputs.get('INPUT_INFO'))
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
                                 'INPUTS_DICT' : inputs
                                 })


    # RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
    nworker = inputs['NWORKER']
    if nworker == -1:
        nworker = cpu_count()

    if nworker == 0:
        for dict_ in mega_list:
            run_fssh(dict_)
    else:
        pool = Pool(nworker)
        pool.map( run_fssh, mega_list)









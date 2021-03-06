#!/usr/bin/python
# standard module
import numpy
from itertools import product as iterprod


# custom modules
from utils_task import *


# SET-UP THE SYSTEM
paths, nworker, archer = set_up(sys.argv)


info = {
    #################### CAN BE CHANGED ###############################################
    #'FILE_INIT': ['run-sample-%s' % x for x in range(1)],  # NAME OF THE RUN OF INITIALIZATION
    'FILE_INIT': ['run-sample-%s' % x for x in range(1)],  # NAME OF THE RUN OF INITIALIZATION
    'NUMBER_CONFIG': 1,
    'NUMBER_REPEAT': 1,
    'LENGTH_FS':   10,  # LENGTH IN FS
    'PRINTING_FREQUENCY_FAST' : 1, # Print every N fs
    'AOM_RADIUS' : 3.0,
  #  'NAME': 'anthracene-fssh',
  #  'SEED': 717063212,
    ###################################################################################
    'INITIALIZATION': 'DIABATIC',
    'LIGHT': False,
}
info['PRINTING_FREQUENCY_SLOW'] = info['LENGTH_FS']
cp2k_param = [
    #################### CAN BE CHANGED ###############################################
    ['PROPAGATION', 'FSSH'],  # METHOD OF PROPAGATION: FSSH OR BORN_OPPENEHIMER
    ['SCALING', 0.06685],  # SCALING FACTOR IN HARTREE (C = 1.819 eV)
    ['TIMESTEP', 0.5],  # TIMESTEP IN FS
    ['REPEAT'] + range(info.get('NUMBER_REPEAT')),  # NUMBER OF FSSH RUN PER STARTING POINT
    #['FILE_INIT'] + info['FILE_INIT'],
    ###################################################################################
    ['DECOHERENCE_CORRECTIONS', 'DAMPING'],
    ['SPURIOUS_TRANSFER_CORR','T'],   
    ['REORDERING_STATES_USING_OVERLAP', 'T'], 
    ['DECO_TIME', 'FORCES_BASED'],
    ['THRESHOLD_TAU_FORCES',  1E-20],
    ['TEMPERATURE_FG_WIDTH',  298],
    ['METHOD_RESCALING', 'NACV'],
    ['METHOD_ADIAB_NACV', 'FAST'],
    ['METHOD_REVERSAL', 'ALWAYS'],
    ['EDC_E0', 0.1],
    ['ELECTRONIC_STEPS', 5],
    ['TEMPLATE_FILE', 'FSSH_PBC_CRYSTAL.template'],
    ['FORCEFIELD_FILE', 'ANTHRACENE_HOLE_FF.prm'],
    ['AOM_COEFF', 'ANTHRACENE_HOLE_AOM.inc'],
    ['INITIALIZATION', 'DIABATIC'],
    ['PRINT_MORE', 'F'],
    ['RUNLOG', 'OFF']
]

# BUILD THE MEGA_LISTS
mega_list = []
for init in info['FILE_INIT']:
    system_info = (InputFile('%s/system.info' % init).dict)
    system_info.update({'AOM_RADIUS' : info['AOM_RADIUS']})
    system_info.update({'FILE_INIT' : init})
    system_info = add_list_activated(system_info, init)
    list_config_init = range(0, system_info['NPROD_INIT'], system_info['NPROD_INIT'] / system_info['NCONFIG_INIT'])
    cp2k_param_here = list(cp2k_param)
    cp2k_param_here.append(['INIT'] + [list_config_init[x] for x in
     range(0, len(list_config_init), len(list_config_init) / info['NUMBER_CONFIG'])])
    second_list = [sublist[1:] for sublist in cp2k_param_here]
    total_list = list(iterprod(*second_list))
    for sublist in total_list:
        subdict = {}
        for index in range(len(sublist)):
            subdict.update({
                cp2k_param_here[index][0]: sublist[index]
            })
        subdict.update(system_info)
        subdict.update({'ARCHER' : archer})
        mega_list.append(subdict)
for ndir in range(len(mega_list)):
    mega_list[ndir].update({'NDIR': ndir,
                            'PATHS_DICT': paths,
                            'INPUTS_DICT': info
                            })



# RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
if nworker < 2:
    for cp2k_info in mega_list:
        run_fssh(cp2k_info)
elif nworker > 1:
    from multiprocessing import Pool
    pool = Pool(nworker)
    pool.map(run_fssh, mega_list)



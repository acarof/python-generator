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
    'FILE_INIT_TITLE': 'TASK345-DIMER-BO-Long200ps_0.1fs_0.2lambda-171002-3ea76b49287023fbfed1817018642dd4',  # NAME OF THE RUN OF INITIALIZATION
    'NUMBER_CONFIG': 100,
    'NUMBER_REPEAT': 10,
    'LENGTH_FS':   10,  # LENGTH IN FS
    'AOM_RADIUS' : 3.0,
    'SCALING' : [0.03, 0.02, 0.01, 0.008, 0.003, 0.0005],
    #'SCALING' : [0.03, 0.02],
    'STATE' :   [1, 2],
    'NPROD_INIT' : 2000000,
    'NCONFIG_INIT' : 1000,
    ###################################################################################
    'INITIALIZATION': 'ADIABATIC',
    'LIGHT': False,
}
cp2k_param = [
    #################### CAN BE CHANGED ###############################################
    ['PROPAGATION', 'FSSH'],  # METHOD OF PROPAGATION: FSSH OR BORN_OPPENEHIMER
    #['SCALING'] + info['SCALING'],
    ['TIMESTEP', 0.5],  # TIMESTEP IN FS
    #['REPEAT'] + range(info.get('NUMBER_REPEAT')),  # NUMBER OF FSSH RUN PER STARTING POINT
    #['FILE_INIT'] + info['FILE_INIT'],
    ###################################################################################
    ['DECOHERENCE_CORRECTIONS', 'DAMPING'],
    ['DECO_TIME', 'FORCES_BASED'],
    ['THRESHOLD_TAU_FORCES',  1E-20],
    ['TEMPERATURE_FG_WIDTH',  298],
    ['METHOD_RESCALING', 'NACV'],
    ['METHOD_ADIAB_NACV', 'FAST'],
    ['METHOD_REVERSAL', 'ALWAYS'],
    ['EDC_E0', 0.1],
    ['ELECTRONIC_STEPS', 5],
    ['TEMPLATE_FILE', 'FSSH_PBC_CRYSTAL.template'],
    ['FORCEFIELD_FILE', 'ANTRACENE_FF.prm'],
    ['INITIALIZATION', 'DIABATIC'],
]


number_traj = info['NUMBER_CONFIG']*info['NUMBER_REPEAT']

#DICTIONARY FOR EQUILIBRIUM POPULATION
dict_for_equilibrium = {}
pop = {}
enough = {}
with open('initial/Simulation-1.dat' ) as file_dict:
    for line in file_dict.readlines():
        if '#' in line:
            pass
        else:
            dict_for_equilibrium[
                float(line.split()[0])
            ] = [float(x) for x in line.split()[1:]]
print dict_for_equilibrium
for scaling in info['SCALING']:
    pop[scaling] = [0] * len(info['STATE'])
    enough[scaling] = [True] * len(info['STATE'])
    pop[scaling][0] = number_traj
    for state in info['STATE']:
        if state != 1:
            round_ = round(dict_for_equilibrium[scaling][state - 1] * number_traj)
            if round_ >= 1:
                pop[scaling][state  - 1] = round_
            else:
                pop[scaling][state - 1] = 1.0
                enough[scaling][state - 1] = False
            pop[scaling][0] += - pop[scaling][state - 1]
print pop
print enough

# BUILD THE MEGA_LISTS
list_config_init = range(0, info['NPROD_INIT'], info['NPROD_INIT'] / info['NCONFIG_INIT'])
print "List config init: ", len(list_config_init)
mega_list = []
for scaling in info['SCALING']:
    print "Scaling: ", scaling
    for state in info['STATE']:
        init = 'initial/state-%s-scaling-%s-%s' % (state, scaling, info['FILE_INIT_TITLE'])
        system_info = (InputFile('%s/system.info' % init).dict)
        system_info['FILE_INIT'] = init
        system_info.update({'AOM_RADIUS' : info['AOM_RADIUS']})
        system_info.update({'NMOL_UNIT' : 1})
        system_info.update({'SYSTEM' : 'OS_SOLVENT'})
        system_info['COORD_CHARGE'] = [0, 1, 1]
        system_info['ELECTROSTATICS'] = False
        system_info['FORCEFIELD_FILE'] = 'ETHYLENE_NE_FF.inc'
        for letter in ['A', 'B', 'C']:
            system_info['VECT_%s' % letter] = system_info['VECT%s' % letter]
        system_info = add_list_activated(system_info, init)
        system_info['FIRST_ADIABAT'] = state
        print pop[scaling][state - 1]
        number_config = int(round( pop[scaling][state - 1] / info['NUMBER_REPEAT'])) + 1
        print "number config: ", number_config
        #step = info['NPROD_INIT'] // number_config
        #final = info['NPROD_INIT']
        cp2k_param_here = list(cp2k_param)
        use_configs = [list_config_init[x] for x in
                range(0, len(list_config_init), len(list_config_init) / number_config ) ]
        if number_config != 1:
            use_configs = use_configs[0:number_config-1]
        print "use configs", len(use_configs)
        cp2k_param_here.append(['INIT'] + use_configs )

        if pop[scaling][state - 1]  > info['NUMBER_REPEAT']:
            cp2k_param_here.append( ['REPEAT'] + range(info.get('NUMBER_REPEAT')) )
        else:
            cp2k_param_here.append(['REPEAT'] + [0] * int(pop[scaling][state - 1]) )
        cp2k_param_here.append( ['SCALING', scaling ])
        second_list = [sublist[1:] for sublist in cp2k_param_here]
        total_list = list(iterprod(*second_list))
        print "length total list", len(total_list)
        for sublist in total_list:
            subdict = {}
            for index in range(len(sublist)):
                subdict.update({
                    cp2k_param_here[index][0]: sublist[index]
                })
            subdict.update(system_info)
            subdict.update({'ARCHER' : archer})
            mega_list.append(subdict)

print "length mega_list", len(mega_list)
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



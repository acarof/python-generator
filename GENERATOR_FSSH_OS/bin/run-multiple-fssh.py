#!/usr/bin/python
# standard module
import numpy
from itertools import product as iterprod


# custom modules
from utils_task import *
import stat

# SET-UP THE SYSTEM
paths, nworker, archer = set_up(sys.argv)


# muliple_info is a list of dictionnary
# Each dictionary will generate one or several FSSH runs with the parameters defined within
multiple_info = [{}]  # Default dictionnary

multiple_info = [
    {'FILE_INIT' : 'run-sample-0',
     'TIMESTEP'  : 0.1,
     'NEW_DIR'   : 'CHEDDAR_MAN'
     },
#    {'FILE_INIT': 'run-sample-0',
#     'TIMESTEP': 0.5,
#    'NEW_DIR': 'NEW_TEST2'
#},
    #{'FILE_INIT': 'run-sample-0',
    # 'SCALING': 0.11111},
]

info = {
    ###################################################################################
    #### TRAJECTORIES INFORMATION ####
    # 'FILE_INIT' : 'run-sample-0',                    # path to the run where initial config can be found
    #'TIMESTEP': 0.5,                                   # Timestep (fs)
    'NUMBER_CONFIG': 1,                                # number of initial config
    'NUMBER_REPEAT': 1,                                # number of repetition for each config
    'LENGTH_FS': 10,  # LENGTH IN FS                   # length of the run
    'PRINTING_FREQUENCY_FAST': 1,                      # Print every N fs for FSSH properties (coefficient, hamitlonian)

    ###################################################################################
    #### SYSTEM INFORMATION ####
    'FORCEFIELD_FILE': 'ANTHRACENE_HOLE_HUI_FF.prm',       # Force field file
    'FORCEFIELD_LJ_ONLY': 'ANTHRACENE_ONLYLJ_FF.prm',  # Force field file with only LJ parameters
    'AOM_COEFF': 'ANTHRACENE_HOLE_HUI_AOM.inc',            # AOM P-PI coefficient
    'TOPOLOGY_NOBOND' : 'ANTHRACENE_NOBOND_TOPOLOGY.inc',    # File with molecule topology without bond
    'SCALING': 0.1136,                                # Scaling factor for AOM method (C = 1.819 eV)
    'AOM_RADIUS': 3.0,                                 # Radius to define the chain of molecule part of FOBSH

    ###################################################################################
    #### RUN INFORMATION ####
    'ARCHER': archer,                                  # Logical to know if we have to run the simulation or not
    #  'NAME': 'anthracene-fssh',                      # Name of the runs
    #  'SEED': 717063212,                              # Seed of the runs, if not present, a random value is chosen
    'LIGHT': False,                                    # Remove some files for production run
    'PRINT_MORE' : 'T',                                # Level of printing for FSSH information ( T or F)
    'RUNLOG'     : 'LOW',                              # Level of printing for MD information (OFF or LOW)
    'DO_SPEEDUP_LJ' : 'T',                             # Use LJ speed-up (T or F)
    'DO_SPEEDUP_INTRA' : 'T',                          # Use intra speed-up (T or F)

    ###################################################################################
    #### SURFACE HOPPING INFORMATION ####
    'PROPAGATION': 'FSSH',                             # Method of propagation: FSSH or BORN_OPPENEHIMER
    'INITIALIZATION': 'DIABATIC',                      # What kind of initialization: ADIABATIC or DIABATIC
    'ELECTRONIC_STEPS': 5,                             # How many electronic steps in one nuclear timestep (for electronic propagation)
    'TEMPLATE_FILE': 'FSSH_PBC_CRYSTAL.template',      # Template for FSSH simulation

    ###################################################################################
    #### DECOHERENCE INFORMATION ####
    'DECOHERENCE_CORRECTIONS': 'DAMPING',              # Decoherence correction: DAMPING, NO_DECO_CORR, etc.
    'SPURIOUS_TRANSFER_CORR' : 'T',                    # Apply spurious transfer correction : True ('T') or False ('F')
    'REORDERING_STATES_USING_OVERLAP': 'T',            # Apply reordering of the state with overlap : True ('T') or False ('F')
    'DECO_TIME': 'FORCES_BASED',                       # Which deco time for decoherence: FORCES_BASED, EDC, etc.
    'THRESHOLD_TAU_FORCES': 1E-20,                     # If DECO_TIME is FORCES_BASED, we have a threshold to speed up the calculation
  #  'TEMPERATURE_FG_WIDTH': 298,                       # If DECO_TIME is FORCES_BASED, we need the temperature as parameters
    'EDC_E0': 0.1,                                     # If DECO_TIME is EDC, parameter E0

    ###################################################################################
    #### RESCALING/REVERSAL INFORMATION ####
    'METHOD_RESCALING': 'NACV',                        # Rescaling of velocities after a successful hop: NACV, SIMPLE, QSYS_SIMPLE
    'METHOD_ADIAB_NACV': 'FAST',                       # If rescaling is NACV, method to calculate them: FAST, TOTAL
    'METHOD_REVERSAL': 'ALWAYS',                       # Reversal of velocities after a rejected hop: NEVER, ALWAYS, TRHULAR, SUBOTNIK



}
info['PRINTING_FREQUENCY_SLOW'] = info['LENGTH_FS']


def my_update(dict1, dict2):
    if (len(set(dict1.keys()) & set(dict2.keys())) != 0):
        print "At least one key is common to info and multiple info:"
        print set(dict2.keys()) & set(dict1.keys())
        raise SystemExit
    else:
        dict1.update(info)


# BUILD THE MEGA_LISTS
mega_list = []
for key_old in multiple_info:
    key = dict(key_old)
    my_update(key, info)
    init = key['FILE_INIT']

    system_info = (InputFile('%s/system.info' % init).dict)
    system_info.update({'AOM_RADIUS' : key['AOM_RADIUS']})
    system_info.update({'FILE_INIT' : init})
    system_info = add_list_activated(system_info, init)
    list_config_init = range(0, system_info['NPROD_INIT'], system_info['NPROD_INIT'] / system_info['NCONFIG_INIT'])
    list_init = [list_config_init[x] for x in
     range(0, len(list_config_init), len(list_config_init) /key['NUMBER_CONFIG'])]
    list_repeat = range(key['NUMBER_REPEAT'])


    for repeat in list_repeat:
        for number in list_init:
            new_key = dict(key)
            system_info.update({ 'INIT' : number})
            new_key.update(system_info)
            #print new_key
            mega_list.append(new_key)


move_list = []
for ndir in range(len(mega_list)):
    mega_list[ndir].update({'NDIR': ndir,
                            'PATHS_DICT': paths,
                            'INPUTS_DICT': info
                            })
    move_list.append( mega_list[ndir].get('NEW_DIR', '.' ))


# RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
if nworker < 2:
    for cp2k_info in mega_list:
        run_fssh(cp2k_info)
elif nworker > 1:
    from multiprocessing import Pool
    pool = Pool(nworker)
    pool.map(run_fssh, mega_list)


target_ndir = 0
target_list = []
for ndir in range(len(mega_list)):
    target = move_list[ndir]
    if target not in target_list:
        target_list.append(target)
    if target != '.':
        if not (os.path.isdir(target)):
            os.mkdir(target)
        # WARNING: we assume name is run-fssh-
        os.system( 'mv run-fssh-%s %s/run-fssh-%s' % (ndir, target, target_ndir))
    target_ndir += 1
    if target_ndir >= move_list.count(target):
        target_ndir = 0


def create_x(results, name):
    with open(name, 'w') as file_:
        file_.write(results)
    st = os.stat(name)
    os.chmod(name, st.st_mode | stat.S_IEXEC)


results = """
#!/bin/bash

list_run=(%s)

for dir in "${list_run[@]}";
do
    cp job_cp2k_taskfarm_fssh.pbs ${dir}/
    cp cpuid.x ${dir}/
    cp task-for-fssh-taskfarm.sh ${dir}/
    cp -r topologies ${dir}/
    cd ${dir}
    pwd
    qsub job_cp2k_taskfarm_fssh.pbs &
    cd ..
done
""" % "  ".join(["\"%s\"" % target for target in target_list])
create_x(results, name="do_multiple_fssh.sh")



results = """
#!/bin/bash


list_run=(%s)

for dir in "${list_run[@]}";
do
    cp job_analyse_general.pbs ${dir}/
    cp -r scripts/ ${dir}/
    cd ${dir}
    pwd
    qsub job_analyse_general.pbs &
    cd ..
done
""" % "  ".join(["\"%s\"" % target for target in target_list])
create_x(results, name="do_multiple_analysis.sh")


results = """
#!/bin/bash


list_run=(%s)

for dir in "${list_run[@]}";
do
    cp job_archive_bash.pbs ${dir}/
    cp -r scripts/ ${dir}/
    cd ${dir}
    pwd
    qsub job_archive_bash.pbs &
    cd ..
done
""" % "  ".join(["\"%s\"" % target for target in target_list])
create_x(results, name="do_multiple_archive.sh")



results = """
#!/bin/bash


list_run=(%s)

for dir in "${list_run[@]}";
do
    cp job_cut_bash.pbs ${dir}/
    cp -r scripts/ ${dir}/
    cd ${dir}
    pwd
    qsub job_cut_bash.pbs &
    cd ..
done
""" % "  ".join(["\"%s\"" % target for target in target_list])
create_x(results, name="do_multiple_cut.sh")

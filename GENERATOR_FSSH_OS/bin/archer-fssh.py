#!/usr/bin/python
# The jobname
#PBS -N generator

# The total number of nodes for your job.
#PBS -l select=1
nworker = 0

# The walltime needed
#PBS -l walltime=1:0:0

# Set the budget to charge the job to (change "budget" to your code)
#PBS -A e358

# standard module
import numpy
from itertools import product as iterprod


# custom modules
from utils_task import *
from find_crystal_bb import find_molecules


# FIND CP2K PATH
paths = {}
def create_cp2k_path():
    print "The file cp2k.path doesn't exist, please provide the path for CP2K:"
    path = raw_input('> ')
    cp2k_path = open('bin/cp2k.path', 'w')
    cp2k_path.write(path)
    cp2k_path.close()
    return open('bin/cp2k.path', 'r')
try:
    cp2k_path = open('bin/cp2k.path', 'r')
except:
    cp2k_path = create_cp2k_path()
paths.update(
    {'cp2k' : cp2k_path.readline()})


task_info = {
    #################### CAN BE CHANGED ###############################################
    'KIND_RUN': 'TONAME',  # NAME OF YOUR RUN
    'FILE_INIT': 'initial/from-GENERATOR_FSSH_OS-temp-1',  # NAME OF THE RUN OF INITIALIZATION
    'NUMBER_CONFIG': 1,
    'NUMBER_REPEAT': 1,
    'LENGTH_FS': 1,  # LENGTH IN FS
    ###################################################################################
    'INITIALIZATION': 'DIABATIC',
    'LIGHT': True,
}
system_info = (InputFile('%s/system.info' % task_info['FILE_INIT']).dict)
system_info.update({
    'AOM_RADIUS': 3.0
})
list_config_init = range(0, system_info['NPROD_INIT'], system_info['NPROD_INIT'] / system_info['NCONFIG_INIT'])
cp2k_param = [
    #################### CAN BE CHANGED ###############################################
    ['PROPAGATION', 'FSSH'],  # METHOD OF PROPAGATION: FSSH OR BORN_OPPENEHIMER
    ['SCALING', 0.06685],  # SCALING FACTOR IN HARTREE (C = 1.819 eV)
    ['TIMESTEP', 0.5],  # TIMESTEP IN FS
    ['REPEAT'] + range(task_info.get('NUMBER_REPEAT')),  # NUMBER OF FSSH RUN PER STARTING POINT
    ###################################################################################
    ['INIT'] + \
    [list_config_init[x] for x in
     range(0, len(list_config_init), len(list_config_init) / task_info['NUMBER_CONFIG'])],  # STARTING POINTS
    ['DECOHERENCE_CORRECTIONS', 'DAMPING'],
    ['DECO_TIME', 'FORCES_BASED'],
    ['THRESHOLD_TAU_FORCES',  2.0],
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

# FIND ACTIVE MOLECULES WITH GUIDO'S TOOL
length = 0.0
for x, y in zip(system_info['ABC'], system_info['DIRECTION']):
    length += (x * y) ** 2
length = numpy.sqrt(length) * (system_info['NUMBER_MOL_ACTIVE'] - 1)
system_info.update({
    'LIST_ACTIVATED': find_molecules(
        coord_first=system_info['COORD_FIRST'],
        size_crystal=system_info['SIZE_CRYSTAL'],
        length=length,
        vector=system_info['DIRECTION'],
        radius_aom=system_info['AOM_RADIUS'],
        psf_file='%s/input-1.psf'  % task_info['FILE_INIT'],
        xyz_file='%s/crystal.xyz'  % task_info['FILE_INIT'],
        nmol_unit=system_info['NMOL_UNIT']

    ),
    'FIRST_DIABAT': int(numpy.ceil(system_info['NUMBER_MOL_ACTIVE'] / 2))
})

print system_info['LIST_ACTIVATED']
print "First diabat: ", system_info['FIRST_DIABAT']
# system_info['LIST_ACTIVATED'] = [1, 2, 3]
# exit()

# BUILD THE MEGA_LISTS
second_list = [sublist[1:] for sublist in cp2k_param]
total_list = list(iterprod(*second_list))
mega_list = []
for sublist in total_list:
    subdict = {}
    for index in range(len(sublist)):
        subdict.update({
            cp2k_param[index][0]: sublist[index]
        })
    mega_list.append(subdict)

# SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
paths.update({'bucket': os.getcwd()})
# PATHS CONTAINS ALL THE PATHS
for directory in ['bin', 'initial', 'scripts', 'structures', 'tasks', 'templates', 'tools', 'topologies']:
    dir = Dir(directory, paths)
    dir.checkdir()

# PREPARE THE MEGA_LIST FOR POOL
for ndir in range(len(mega_list)):
    mega_list[ndir].update({'NDIR': ndir,
                            'PATHS_DICT': paths,
                            'INPUTS_DICT': task_info
                            })
    mega_list[ndir].update(system_info)

# RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
if nworker == 0:
    for cp2k_info in mega_list:
        run_fssh(cp2k_info)
else:
    from multiprocessing import Pool
    pool = Pool(nworker)
    pool.map(run_fssh, mega_list)

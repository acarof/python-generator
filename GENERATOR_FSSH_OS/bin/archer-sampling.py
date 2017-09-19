#!/usr/bin/python
# The jobname
#PBS -N generator

# The total number of nodes for your job.
#PBS -l select=1
nworker = 1

# The walltime needed
#PBS -l walltime=1:0:0

# Set the budget to charge the job to (change "budget" to your code)
#PBS -A e358

# standard modules
import numpy

# custom modules
from utils_task import *
from find_crystal_bb import find_crystal_bb


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
    {'cp2k' : cp2k_path.readline()}
)

def prepare_system_info(dict, path):
    with open('%s/system.info' % path, 'w') as file_:
        for key in dict:
            if key not in ['TEMPLATE_FILE', 'FILE_CRYSTAL', 'FILE_UNIT', 'TIMESTEP']:
                if isinstance(dict[key], (list, tuple)):
                    file_.write('%s    %s\n' % (key, '  '.join(map(str, dict[key]))))
                else:
                    file_.write('%s    %s\n' % (key, dict[key]))



task_info = {
    #################### CAN BE CHANGED ###############################################
    'KIND_RUN': 'TONAME',   # NAME OF YOUR RUN
    'NEQ': 1,              # NUMBER OF TIMESTEP FOR EQUILIBRATION (NVT)
    'NPROD': 1,            # NUMBER OF TIMESTEP FOR PRODUCTION (NVE),
    'TEMPERATURE' : [1],
    'NCONFIG': 1,            # NUMBER OF PRINTED SNAPSHOT
    'ARCHER' : False,
    'NCORE'    : 1
    ##################################################################################
}
system_info = {
    #################### CAN BE CHANGED ###############################################
    'NUMBER_MOL_ACTIVE': 2,                 # NUMBER OF ACTIVE MOLECULES
    'DIRECTION': [1, 0, 0],                 # DIRECTION TO PROPAGATE THE CHARGE
    'RCUT': 8 ,                              # VDW RCUT
    ###################################################################################
    'SYSTEM': 'PBC_CRYSTAL',                                # (do not change)
    'MOL_NAME'            : 'ANTRACENE'          ,          # NAME OF THE MOLECULE
    'FILE_UNIT'           : 'ant_unitcell.xyz'     ,        # NAME OF THE .xyz FILE WITH THE UNITCELL
    'FILE_CRYSTAL'        : 'crystal.xyz',                  # NAME OF THE .xyz FILE TO PRINT THE CRYSTAL
    'ABC'                 : [8.562, 6.038, 11.184],         # ABC OF THE UNITCELL
    'ALPHA_BETA_GAMMA'    : [90.0, 124.70, 90.0],           # ALPHA_BETA_GAMMA OF THE UNITCELL
    'STARTING_POINT'      : [0.0,   0.0,    0.0],           # (do not change)
    'NATOM_MOL'           : 24,                             # NUMBER OF ATOMS PER MOLECULES
    'NMOL_UNIT'           : 2,                              # NUMBER OF MOLECULES PER UNIT_CELL
}
cp2k_info = {
    #################### CAN BE CHANGED ###############################################
     'TIMESTEP'      :       0.5,                           # TIMESTEP IN FS
    ###################################################################################
     'ALPHA'         :    3.5 / system_info['RCUT'],        # ALPHA FOR EWALD
     'TEMPLATE_FILE' : 'FIST_PBC_CRYSTAL.template',         # (do not change)
     'FORCEFIELD_FILE': 'ANTRACENE_FF.prm',                 # FORCEFIELD
}


# FIND CRYSTAL SIZE WITH GUIDO'S TOOL
length = 0.0
for x,y in zip(system_info['ABC'], system_info['DIRECTION']):
    length += (x*y)**2
length= numpy.sqrt(length)*system_info['NUMBER_MOL_ACTIVE']
size_crystal, coord_first = find_crystal_bb(
        abc=system_info['ABC'] + system_info['ALPHA_BETA_GAMMA'],
        start= system_info['STARTING_POINT'],
        length = length,
        vector = system_info['DIRECTION'],
        radius = system_info['RCUT']
    )
print "First molecules of the chain: ", coord_first
coord_charge = [ int(x) for x in (numpy.array(coord_first) \
               + (numpy.ceil(system_info['NUMBER_MOL_ACTIVE']/2))* numpy.array(system_info['DIRECTION']))]
print "Coord charge: ", coord_charge
system_info.update({
    'SIZE_CRYSTAL' : size_crystal,
    'COORD_FIRST' : coord_first,
    'COORD_CHARGE' : coord_charge,
    'LENGTH'       : length
})
cp2k_info.update({
    'GMAX' : [ int(x*y) + 1 for (x,y) in zip(
        system_info['ABC'], system_info['SIZE_CRYSTAL']) ]
})
print 'GMAX = ', cp2k_info['GMAX']



# SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
paths.update({'bucket': os.getcwd()})
# PATHS CONTAINS ALL THE PATHS
for directory in ['bin', 'initial', 'scripts', 'structures', 'tasks', 'templates', 'tools', 'topologies']:
    dir = Dir(directory, paths)
    dir.checkdir()


ndir = 0
generate_initial_structure(system_info, paths)
mega_list = []

for temperature in task_info['TEMPERATURE']:
    output = Dir('initial/from-%s-temp-%s' % (paths['bucket'].split('/')[-1], temperature), paths = paths, target = 'output_here-temp-%s' % temperature)
    output.rm_mkdir()
    mega_list.append(
        {'TEMPERATURE' : temperature,
         'NDIR'        : ndir,
         'SYSTEM_INFO'  : system_info,
         'CP2K_INFO' : cp2k_info,
         'PATHS' : paths,
         'TASK_INFO' : task_info
          }
    )
    ndir += 2


def do_run(dict_):
    temperature = dict_['TEMPERATURE']
    ndir = dict_['NDIR']
    system_info = dict_['SYSTEM_INFO']
    cp2k_info = dict_['CP2K_INFO']
    paths     = dict_['PATHS']
    task_info = dict_['TASK_INFO']


    previous_dir = run_fist(system_info, cp2k_info, paths, steps=task_info['NEQ'],
                                  ndir=ndir, restart_info=None, velocities=False, ensemble='NVT',
                                  TEMPERATURE=temperature, nconfig=task_info['NCONFIG'],
                                  archer=task_info['ARCHER'], nworker=max(task_info['NCORE'], 1), name = 'eq')
    print previous_dir
    ndir += 1

    restart_info = {
        'RESTART_DIR': previous_dir,
        'CONFIG': task_info['NEQ']
    }
    previous_dir = run_fist(system_info, cp2k_info, paths, steps=task_info['NPROD'],
                                  ndir=ndir, restart_info=restart_info, velocities=True, ensemble='NVE',
                                  TEMPERATURE=temperature, nconfig=task_info['NCONFIG'],
                                  archer=task_info['ARCHER'], nworker=max(task_info['NCORE'], 1), name = 'sample')

    os.system('cp %s/run-pos-1.xyz %s' % (previous_dir, paths['output_here-temp-%s' % temperature]))
    os.system('cp %s/run-vel-1.xyz %s' % (previous_dir, paths['output_here-temp-%s' % temperature]))
    os.system('cp %s/input-1.psf %s' % (previous_dir, paths['output_here-temp-%s' % temperature]))
    os.system('cp %s/crystal.xyz %s' % (paths['crystal'], paths['output_here-temp-%s' % temperature]))
    system_info.update({'TEMPERATURE': temperature})
    system_info.update({'NPROD_INIT' : task_info['NPROD']})
    system_info.update({'NCONFIG_INIT' : task_info['NCONFIG']})
    prepare_system_info(system_info, paths['output_here-temp-%s' % temperature])


# RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
if nworker == 0:
    for cp2k_info in mega_list:
        do_run(cp2k_info)
else:
    from multiprocessing import Pool
    pool = Pool(nworker)
    pool.map( do_run, mega_list)

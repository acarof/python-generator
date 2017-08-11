#!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy

# custom modules
from utils import Dir, Bucket
from utils_task import *
from find_crystal_bb import find_crystal_bb


def prepare_system_info(dict, path):
    with open('%s/system.info' % path, 'w') as file_:
        for key in dict:
            if key not in ['TEMPLATE_FILE', 'FILE_CRYSTAL', 'FILE_UNIT', 'TIMESTEP']:
                if isinstance(dict[key], (list, tuple)):
                    file_.write('%s    %s\n' % (key, '  '.join(map(str, dict[key]))))
                else:
                    file_.write('%s    %s\n' % (key, dict[key]))

def prepare_future_task(path):
    with open('%s/tasks/task-fssh-pbc-crystal/task.py' % path, 'w') as fileout, open('tasks/task-fssh-pbc-crystal/task.py') as filein:
        for line in filein.readlines():
            if "'FILE_INIT':" in line:
                result = "        'FILE_INIT': '%s',           # NAME OF THE RUN OF INITIALIZATION\n" % path.split('/')[-2][5:]
                fileout.write(result)
            else:
                fileout.write(line)


def main(inputs, paths):
    print """
    START THE DIABAT OF SAMPLING OF AN ORGANIC CRYSTAL
    """

    task_info = {
        #################### CAN BE CHANGED ###############################################
        'KIND_RUN': 'TONAME',   # NAME OF YOUR RUN
        'NEQ': 20,              # NUMBER OF TIMESTEP FOR EQUILIBRATION (NVT)
        'NPROD': 100,            # NUMBER OF TIMESTEP FOR PRODUCTION (NVE),
        'TEMPERATURE' : [100],
        'NCONFIG': 10,            # NUMBER OF PRINTED SNAPSHOT
        'PARALLEL' : True
#        'TEMPERATURE' : [100, 140, 180, 220, 260, 300 ]
        ##################################################################################
    }


    system_info = {
        #################### CAN BE CHANGED ###############################################
        'NUMBER_MOL_ACTIVE': 12,                 # NUMBER OF ACTIVE MOLECULES
        'DIRECTION': [0, 1, 0],                 # DIRECTION TO PROPAGATE THE CHARGE
        'RCUT': 8 ,                              # VDW RCUT
        ###################################################################################
        'SYSTEM': 'NEUTRAL_CRYSTAL',                                # (do not change)
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
         'TEMPLATE_FILE' : 'FIST_PBC_CRYSTAL.template',         # (do not change)
         'FORCEFIELD_FILE': 'ANTRACENE_FF.prm',                 # FORCEFIELD
    }


    # FIND CRYSTAL SIZE WITH GUIDO'S TOOL
    length = 0.0
    for x,y in zip(system_info['ABC'], system_info['DIRECTION']):
        length += (x*y)**2
    length= numpy.sqrt(length)*system_info['NUMBER_MOL_ACTIVE']
    size_crystal, coord_charge = find_crystal_bb(
            abc=system_info['ABC'] + system_info['ALPHA_BETA_GAMMA'],
            start= system_info['STARTING_POINT'],
            length = length,
            vector = system_info['DIRECTION'],
            radius = system_info['RCUT']
        )
    system_info.update({
        'SIZE_CRYSTAL' : size_crystal,
        'COORD_CHARGE' : coord_charge,
        'LENGTH'       : length
    })


    # SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
    inputs.update(task_info)
    bucket = Bucket(inputs)
    bucket.name()
    paths.update({'bucket': bucket.path})
    # PATHS CONTAINS ALL THE PATHS
    for directory in ['bin', 'initial', 'scripts', 'structures', 'tasks', 'templates', 'tools', 'topologies']:
        dir = Dir(directory, paths)
        dir.checkdir()



    output = Dir('output', paths = paths)
    output.rm_mkdir()
    ndir = 0


    for temperature in task_info['TEMPERATURE']:
        output = Dir('output/from-%s-temp-%s' % (paths['bucket'].split('/')[-1], temperature), paths = paths, target = 'output_here')
        output.rm_mkdir()
        output_initial = Dir('output/from-%s-temp-%s/initial/' % (paths['bucket'].split('/')[-1], temperature), paths = paths, target = 'output_initial')
        output_initial.rm_mkdir()
        output_initial = Dir('output/from-%s-temp-%s/initial/from-%s-temp-%s' % (paths['bucket'].split('/')[-1], temperature, paths['bucket'].split('/')[-1], temperature),
                             paths = paths, target = 'output_initial_here')
        output_initial.rm_mkdir()



        generate_initial_structure(system_info, paths)
        ndir, previous_dir = run_fist(system_info, cp2k_info, paths, steps = task_info['NEQ'],
                                      ndir = ndir, restart_info = None, velocities = False, ensemble = 'NVT',
                                      TEMPERATURE=temperature, nconfig = task_info['NCONFIG'],
                                      parallel = task_info['PARALLEL'], nworker = max(inputs['NWORKER'], 1) )
        print ndir, previous_dir

        restart_info = {
            'RESTART_DIR' : 'run-%s' % previous_dir,
            'CONFIG'      : task_info['NEQ']
        }
        ndir, previous_dir = run_fist(system_info, cp2k_info, paths, steps = task_info['NPROD'],
                                      ndir = ndir, restart_info = restart_info, velocities = True, ensemble = 'NVE',
                                      TEMPERATURE=temperature, nconfig = task_info['NCONFIG'],
                                      parallel = task_info['PARALLEL'], nworker = max(inputs['NWORKER'], 1))

        os.system('cp run-%s/run-pos-1.xyz %s' % (previous_dir, paths['output_initial_here']))
        os.system('cp run-%s/run-vel-1.xyz %s' % (previous_dir, paths['output_initial_here']))
        os.system('cp run-%s/input-1.psf %s' % (previous_dir, paths['output_initial_here']))
        os.system('cp %s/crystal.xyz %s' % (paths['crystal'], paths['output_initial_here']))
        system_info.update({ 'TEMPERATURE' : temperature})
        prepare_system_info(system_info, paths['output_initial_here'])
        for directory in ['bin', 'scripts','structures', 'templates', 'tools', 'topologies']:
            os.system('cp -r %s %s' % (directory, paths['output_here']))
        os.system('mkdir -p %s/tasks/task-fssh-pbc-crystal' % paths['output_here'])
        prepare_future_task(paths['output_here'])

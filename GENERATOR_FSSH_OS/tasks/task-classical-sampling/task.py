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

def main(inputs, paths):
    print """
    START THE DIABAT OF SAMPLING OF AN ORGANIC CRYSTAL
    """

    task_info = {
        'KIND_RUN': 'TONAME',
        'NEQ': 20,
        'NPROD': 10,
        'NCONFIG': 3
    }


    system_info = {
        'SYSTEM'              : 'PBC_CRYSTAL'            ,
        'MOL_NAME'            : 'ANTRACENE'          ,
        'FILE_UNIT'           : 'ant_unitcell.xyz'     ,
        #'FILE_UNIT'           : 'ant_c.xyz'     ,
        'FILE_CRYSTAL'        : 'crystal.xyz',
        'ABC'                 : [8.562, 6.038, 11.184],
        'ALPHA_BETA_GAMMA'    : [90.0, 124.70, 90.0],
        #'SIZE_CRYSTAL'        : [4    , 7    ,  2    ] ,
        'DIRECTION'           : [0,     1,      0    ],
        'NUMBER_MOL_ACTIVE'   : 3,
        'STARTING_POINT'      : [0.0,   0.0,    0.0],
        'NATOM_MOL'           : 24,
        'NMOL_UNIT'           : 2,
        'RCUT'                : 8
    }

    cp2k_info = {
        #[ 'TIMESTEP', 0.01, 0.05, 0.1, 0.5],
         'TIMESTEP'      :       0.1,
         'TEMPLATE_FILE' : 'FIST_PBC_CRYSTAL.template',
         'FORCEFIELD_FILE': 'ANTRACENE_FF.prm',
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
    bucket = Bucket(inputs)
    bucket.name()
    paths.update({'bucket': bucket.path})

    output = Dir('output', paths = paths)
    output.rm_mkdir()
    output = Dir('output/from-%s' % paths['bucket'].split('/')[-1], paths = paths, target = 'output_here')
    output.rm_mkdir()
    output_initial = Dir('output/from-%s/initial/' % paths['bucket'].split('/')[-1], paths = paths, target = 'output_initial')
    output_initial.rm_mkdir()
    output_initial = Dir('output/from-%s/initial/from-%s' % (paths['bucket'].split('/')[-1], paths['bucket'].split('/')[-1]),
                         paths = paths, target = 'output_initial_here')
    output_initial.rm_mkdir()



    generate_initial_structure(system_info, paths)
    ndir, previous_dir = run_fist(system_info, cp2k_info, paths, steps = task_info['NEQ'],
                                  ndir = 0, restart_info = None, velocities = False, ensemble = 'NVT')
    print ndir, previous_dir

    restart_info = {
        'RESTART_DIR' : 'run-%s' % previous_dir,
        'CONFIG'      : task_info['NEQ']
    }
    ndir, previous_dir = run_fist(system_info, cp2k_info, paths, steps = task_info['NPROD'],
                                  ndir = ndir, restart_info = restart_info, velocities = True, ensemble = 'NVE')

    os.system('cp run-%s/run-pos-1.xyz %s' % (previous_dir, paths['output_initial_here']))
    os.system('cp run-%s/run-vel-1.xyz %s' % (previous_dir, paths['output_initial_here']))
    os.system('cp run-%s/input-1.psf %s' % (previous_dir, paths['output_initial_here']))
    os.system('cp %s/crystal.xyz %s' % (paths['crystal'], paths['output_initial_here']))
    prepare_system_info(system_info, paths['output_initial_here'])
    for directory in ['bin', 'scripts',  'tasks', 'templates', 'tools', 'topologies']:
        os.system('cp -r %s %s' % (directory, paths['output_here']))
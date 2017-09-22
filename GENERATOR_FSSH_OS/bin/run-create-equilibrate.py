#!/usr/bin/python
# standard modules
import numpy

# custom modules
from utils_task import *
from find_crystal_bb import find_crystal_bb

# SET-UP THE SYSTEM
nworker, archer = find_nworker(sys.argv)
print "nworker is: %s and archer is: %s" % (nworker, archer)
paths = find_cp2k_path()
paths.update({'bucket': os.getcwd()})
for directory in ['bin', 'structures', 'templates', 'topologies']:
    dir = Dir(directory, paths)
    dir.checkdir()




anthracene = {
    #################### CAN BE CHANGED ###############################################
    'NEQ': 10,                   # NUMBER OF TIMESTEP FOR EQUILIBRATION (NVT)
    'NCONFIG' : 1,
    'TEMPERATURE_LIST' : [100, 140, 180, 220, 260, 300],
    'TIMESTEP': 0.5,            # TIMESTEP IN FS

    #################### CAN BE CHANGED ###############################################
    'NUMBER_MOL_ACTIVE': 12,  # NUMBER OF ACTIVE MOLECULES
    'DIRECTION': [0, 1, 0],  # DIRECTION TO PROPAGATE THE CHARGE
    'RCUT': 8,  # VDW RCUT

    ###################################################################################
    'SYSTEM': 'PBC_CRYSTAL',  # (do not change)
    'MOL_NAME': 'ANTRACENE',  # NAME OF THE MOLECULE
    'FILE_UNIT': 'ant_unitcell.xyz',  # NAME OF THE .xyz FILE WITH THE UNITCELL
    'FILE_CRYSTAL': 'crystal.xyz',  # NAME OF THE .xyz FILE TO PRINT THE CRYSTAL
    'ABC': [8.562, 6.038, 11.184],  # ABC OF THE UNITCELL
    'ALPHA_BETA_GAMMA': [90.0, 124.70, 90.0],  # ALPHA_BETA_GAMMA OF THE UNITCELL
    'STARTING_POINT': [0.0, 0.0, 0.0],  # (do not change)
    'NATOM_MOL': 24,  # NUMBER OF ATOMS PER MOLECULES
    'NMOL_UNIT': 2,  # NUMBER OF MOLECULES PER UNIT_CELL
    'FORCEFIELD_FILE': 'ANTRACENE_FF.prm',  # FORCEFIELD
    ##################################################################################
    'TEMPLATE_FILE' : 'FIST_PBC_CRYSTAL.template',         # (do not change)
}
info = anthracene
if info.get('ELECTROSTATICS', True):
    info.update({'ALPHA'         :    3.5 / info['RCUT'] })        # ALPHA FOR EWALD})


if  info['SYSTEM'] == 'PBC_CRYSTAL':
    # FIND CRYSTAL SIZE WITH GUIDO'S TOOL
    length = 0.0
    for x, y in zip(info['ABC'], info['DIRECTION']):
        length += (x * y) ** 2
    length = numpy.sqrt(length) * info['NUMBER_MOL_ACTIVE']
    size_crystal, coord_first = find_crystal_bb(
        abc=info['ABC'] + info['ALPHA_BETA_GAMMA'],
        start=info['STARTING_POINT'],
        length=length,
        vector=info['DIRECTION'],
        radius=info['RCUT']
    )
    print "First molecules of the chain: ", coord_first
    coord_charge = [int(x) for x in (numpy.array(coord_first) \
                                     + (numpy.ceil(info['NUMBER_MOL_ACTIVE'] / 2)) * numpy.array(info['DIRECTION']))]
    print "Coord charge: ", coord_charge
    info.update({
        'SIZE_CRYSTAL': size_crystal,
        'COORD_FIRST': coord_first,
        'COORD_CHARGE': coord_charge,
        'LENGTH': length
    })
    info.update({
        'GMAX': [int(x * y) + 1 for (x, y) in zip(
            info['ABC'], info['SIZE_CRYSTAL'])]
    })
    print 'GMAX = ', info['GMAX']


# GENERATE THE STRUCTURE
output_structure = Dir('structures/from-create-%s' % paths['bucket'].split('/')[-1], paths=paths, target='crystal')
output_structure.rm_mkdir()
generate_initial_structure(info, paths)


# PREPARE MEGA_LIST FOR RUN IN PARALLEL
ndir = 0
mega_list = []
for temperature in info['TEMPERATURE_LIST']:
    mega_list.append(
        {'TEMPERATURE' : temperature,
         'NDIR'        : ndir,
         'INFO'  : info,
         'PATHS' : paths,
         'ARCHER' : archer
          }
    )
    ndir += 1
def do_run(dict_):
    info = dict_['INFO']
    previous_dir = run_fist(dict_['INFO'], dict_['PATHS'], steps=info['NEQ'],
                                  ndir=dict_['NDIR'], restart_info=None, velocities=False, ensemble='NVT',
                                  TEMPERATURE=dict_['TEMPERATURE'], nconfig=info['NCONFIG'],
                                  archer=dict_['ARCHER'], name = 'eq')
    os.system('cp %s/%s %s' % ( paths['crystal'], info['FILE_CRYSTAL'], previous_dir  ))
    info.update({'TEMPERATURE': dict_['TEMPERATURE']})
    info.update({'NEQ' : info['NEQ']})
    prepare_system_info(info, previous_dir)



# RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
if nworker < 2:
    for cp2k_info in mega_list:
        do_run(cp2k_info)
elif nworker > 1:
    from multiprocessing import Pool
    pool = Pool(nworker)
    pool.map( do_run, mega_list)



############################## OTHER INPUT ###################
anthracene = {
    #################### CAN BE CHANGED ###############################################
    'NEQ': 10,                   # NUMBER OF TIMESTEP FOR EQUILIBRATION (NVT)
    'NCONFIG' : 1,
    'TEMPERATURE_LIST' : [100],
    'TIMESTEP': 0.5,            # TIMESTEP IN FS

    #################### CAN BE CHANGED ###############################################
    'NUMBER_MOL_ACTIVE': 12,  # NUMBER OF ACTIVE MOLECULES
    'DIRECTION': [0, 1, 0],  # DIRECTION TO PROPAGATE THE CHARGE
    'RCUT': 8,  # VDW RCUT

    ###################################################################################
    'SYSTEM': 'PBC_CRYSTAL',  # (do not change)
    'MOL_NAME': 'ANTRACENE',  # NAME OF THE MOLECULE
    'FILE_UNIT': 'ant_unitcell.xyz',  # NAME OF THE .xyz FILE WITH THE UNITCELL
    'FILE_CRYSTAL': 'crystal.xyz',  # NAME OF THE .xyz FILE TO PRINT THE CRYSTAL
    'ABC': [8.562, 6.038, 11.184],  # ABC OF THE UNITCELL
    'ALPHA_BETA_GAMMA': [90.0, 124.70, 90.0],  # ALPHA_BETA_GAMMA OF THE UNITCELL
    'STARTING_POINT': [0.0, 0.0, 0.0],  # (do not change)
    'NATOM_MOL': 24,  # NUMBER OF ATOMS PER MOLECULES
    'NMOL_UNIT': 2,  # NUMBER OF MOLECULES PER UNIT_CELL
    'FORCEFIELD_FILE': 'ANTRACENE_FF.prm',  # FORCEFIELD
    ##################################################################################
    'TEMPLATE_FILE' : 'FIST_PBC_CRYSTAL.template',         # (do not change)
}

trimer_ethylene = {
    #################### CAN BE CHANGED ###############################################
    'NEQ': 10,                   # NUMBER OF TIMESTEP FOR EQUILIBRATION (NVT)
    'NCONFIG' : 1,
    'TEMPERATURE_LIST' : [100],
    'TIMESTEP': 0.5,            # TIMESTEP IN FS
    'DENSITY'  : 0.001,
    'SIZE_CRYSTAL' : [3, 1, 1],
    'COORD_CHARGE' : [0, 1, 1], # WARNING, in O-BASIS. FIRST MOLECULE IS [0, 0, 0]
    'SIZE_BOX'     : [30.0, 30.0, 30.0],
    'NUMBER_MOL_ACTIVE': 2,  # NUMBER OF ACTIVE MOLECULES

    ###################################################################################
    'SYSTEM': 'OS_SOLVENT',  # (do not change)
    'MOL_NAME': 'ETHYLENE',  # NAME OF THE MOLECULE
    'SOLVENT' : 'NE',
    'CLOSEST_DIST' : 5,
    'FILE_UNIT': 'ethylene.xyz',  # NAME OF THE .xyz FILE WITH THE UNITCELL
    'FILE_CRYSTAL': 'crystal.xyz',  # NAME OF THE .xyz FILE TO PRINT THE CRYSTAL
    'VECT_A' :  [3.527, 0.784, -0.166],
    'VECT_B' :  [0.000, 0.000,  0.000],
    'VECT_C' :  [0.000, 0.000,  0.000],
    'STARTING_POINT': [0.0, 0.0, 0.0],  # (do not change)
    'NATOM_MOL': 6,  # NUMBER OF ATOMS PER MOLECULES
    'NMOL_UNIT': 1,  # NUMBER OF MOLECULES PER UNIT_CELL
    'FORCEFIELD_FILE': 'ETHYLENE_NE_FF.inc',  # FORCEFIELD
    'ELECTROSTATICS' : False,
    ##################################################################################
    'TEMPLATE_FILE' : 'FIST_PBC_CRYSTAL.template',         # (do not change)
}

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


info = {
    #################### CAN BE CHANGED ###############################################
    'NCONFIG' : 1,             # NUMBER OF TIMESTEP FOR EQUILIBRATION (NVT)
    'NPROD' : 1,
    'FILE_INIT': 'run-eq',  # NAME OF THE RUN OF INITIALIZATION
    'TEMPERATURE_LIST' : [100],
    'TIMESTEP': 0.5,  # TIMESTEP IN FS
    ##################################################################################

    'TEMPLATE_FILE': 'FIST_PBC_CRYSTAL.template',  # (do not change)
}




# PREPARE MEGA_LIST FOR RUN IN PARALLEL
ndir = 0
mega_list = []
for temperature in info['TEMPERATURE_LIST']:
    system_info = (InputFile('%s-%s/system.info' % (info['FILE_INIT'], ndir)).dict)
    info.update(system_info)
    mega_list.append(
        {'TEMPERATURE' : temperature,
         'NDIR'        : ndir,
          'PATHS' : paths,
         'INFO' : info,
         'ARCHER' : archer,
         'RESTART_INFO' :  {
        'RESTART_DIR': 'run-eq-%s' % ndir,
        'CONFIG': info['NEQ']}
          }
    )
    ndir += 1
def do_run(dict_):
    info = dict_['INFO']
    previous_dir = run_fist(system_info=dict_['INFO'],
                                  paths = dict_['PATHS'], steps=info['NPROD'],
                                  ndir=dict_['NDIR'],restart_info=dict_['RESTART_INFO'], velocities=False, ensemble='NVE',
                                  TEMPERATURE=dict_['TEMPERATURE'], nconfig=info['NCONFIG'],
                                  archer=dict_['ARCHER'], name = 'sample')
    os.system('cp %s-%s/crystal.xyz %s' % (info['FILE_INIT'], dict_['NDIR'], previous_dir))
    info.update({'TEMPERATURE': temperature})
    info.update({'NPROD_INIT' : info['NPROD']})
    info.update({'NCONFIG_INIT' : info['NCONFIG']})
    prepare_system_info(info, previous_dir)



# RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
if nworker == 0:
    for cp2k_info in mega_list:
        do_run(cp2k_info)
else:
    from multiprocessing import Pool
    pool = Pool(nworker)
    pool.map( do_run, mega_list)

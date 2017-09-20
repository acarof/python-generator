#!/usr/bin/python
# standard modules
import numpy

# custom modules
from utils_task import *

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
    'NCONFIG' : 100,             # NUMBER OF TIMESTEP FOR EQUILIBRATION (NVT)
    'NPROD' : 100,
    'FILE_INIT': ['run-eq-%s' % x for x in range(1)],  # NAME OF THE RUN OF INITIALIZATION
    'TIMESTEP': 0.5,  # TIMESTEP IN FS

    ##################################################################################
    'TEMPLATE_FILE': 'FIST_PBC_CRYSTAL.template',  # (do not change)
}




# PREPARE MEGA_LIST FOR RUN IN PARALLEL
ndir = 0
mega_list = []
for init in info['FILE_INIT']:
    system_info = (InputFile('%s/system.info' % (init)).dict)
    system_info.update({'TEMPLATE_FILE' : info['TEMPLATE_FILE']})
    mega_list.append(
        {'NDIR'   : ndir,
         'FILE_INIT' : init,
         'PATHS'  : paths,
         'INFO'   : info,
         'SYSTEM' : system_info,
         'ARCHER' : archer,
         'RESTART_INFO' :  {
               'RESTART_DIR': init,
               'CONFIG': system_info['NEQ']}
          }
    )
    ndir += 1
def do_run(dict_):
    info = dict_['INFO']
    system_info = dict_['SYSTEM']
    previous_dir = run_fist(system_info=system_info,
                                  paths = dict_['PATHS'], steps=info['NPROD'],
                                  ndir=dict_['NDIR'],restart_info=dict_['RESTART_INFO'], velocities=True, ensemble='NVE',
                                  TEMPERATURE=system_info['TEMPERATURE'], nconfig=info['NCONFIG'],
                                  archer=dict_['ARCHER'], name = 'sample')
    os.system('cp %s/crystal.xyz %s' % (dict_['FILE_INIT'], previous_dir))
    system_info.update({'NPROD_INIT' : info['NPROD']})
    system_info.update({'NCONFIG_INIT' : info['NCONFIG']})
    prepare_system_info(system_info, previous_dir)



# RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
if nworker == 1:
    for cp2k_info in mega_list:
        do_run(cp2k_info)
elif nworker > 1:
    from multiprocessing import Pool
    pool = Pool(nworker)
    pool.map( do_run, mega_list)

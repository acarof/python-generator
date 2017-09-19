# standard modules
import string, re, struct, sys, math, os, time
import numpy as np
import imp
from collections import deque
from random import seed, randint

# custom modules
from utils import *


def find_nworker(list_):
    try:
        nworker = list_[1]
    except:
        return 0, True
    try:
        nworker = int(nworker)
        return nworker, False
    except:
        if nworker == 'ARCHER':
            return 0, True
        else:
            print "Arguments should be a number or ARCHER"
            raise SystemExit




def prepare_system_info(dict, path):
    with open('%s/system.info' % path, 'w') as file_:
        for key in dict:
            if key not in ['TEMPLATE_FILE', 'FILE_CRYSTAL', 'FILE_UNIT', 'TIMESTEP']:
                if isinstance(dict[key], (list, tuple)):
                    file_.write('%s    %s\n' % (key, '  '.join(map(str, dict[key]))))
                else:
                    file_.write('%s    %s\n' % (key, dict[key]))


def find_cp2k_path():
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
    return {'cp2k' : cp2k_path.readline()}


def generate_initial_structure(system_info, paths):
    system = system_info['SYSTEM']
    if system in ['PBC_CRYSTAL','NEUTRAL_CRYSTAL'] :
        from utils import OSCrystal as Structure
    else:
        print "No structure generator for ", system
        raise SystemExit
    structure = Structure(system_info, paths)
    structure.construct_organic_crystal()



def run_fist(system_info,  paths = {}, steps = 1, ndir = 0, restart_info = None, velocities = True, ensemble = 'NVE',
             TEMPERATURE=300, nconfig = -1, archer= True, nworker = 1, name = None):
    system = system_info['SYSTEM']
    if system == 'PBC_CRYSTAL':
        from utils import FISTOSCrystal as Config
    elif system == 'NEUTRAL_CRYSTAL':
        from utils import FISTOSNeutralCrystal as Config
    else:
        raise SystemExit

    if nconfig == -1:
        my_print = steps
    else:
        my_print = max( steps/nconfig, 1)

    print "GO FOR RUN %d" % ndir
    config = Config(system_info, paths, ENSEMBLE=ensemble, STEPS=steps, RESTART= restart_info, VELOCITIES=velocities,
                    TEMPERATURE=TEMPERATURE, PRINT=my_print, ARCHER = archer, NWORKER = nworker, name = name)
    return config.run(ndir)



def run_fssh( cp2k_info):
    task_info = cp2k_info.get('INPUTS_DICT')
    paths  = cp2k_info.get('PATHS_DICT')

    if task_info['INITIALIZATION'] == 'DIABATIC':
        run_fssh_from_diabat(cp2k_info, task_info, paths)
    else:
        print "NO METHOD IMPLEMENTED FOR INITIALIZATION"
        sys.exit()


def shorten_log_file(ndir):
    with open('run-%s/run.log' % ndir) as oldlog, open('run-%s/new.log' % ndir, 'w') as newlog:
        result = ' '.join([next(oldlog) for x in xrange(100)])
        result += ' '.join(deque(oldlog, 100))
        newlog.write(result)
    os.system('mv run-%s/new.log run-%s/run.log' % (ndir, ndir))
    os.system('rm run-%s/run-r-1.out' % ndir)
    os.system('rm run-%s/run-mix-1.ener' % ndir)
    os.system('rm run-%s/core.*' % ndir)


def run_fssh_from_diabat(cp2k_info, task_info, paths):
    cp2k_info.update({'STEPS' : int(task_info['LENGTH_FS'] / cp2k_info['TIMESTEP'] ) } )
    cp2k_info.update({'PRINT':  int(task_info['LENGTH_FS'] / cp2k_info['TIMESTEP']) } )
    cp2k_info.update({'PRINT_FSSH': int( 1 / cp2k_info['TIMESTEP']) })

    restart_info = {
        'RESTART_DIR' : task_info['FILE_INIT'],
        'CONFIG'      : cp2k_info['INIT']
    }

    seed()
    cp2k_info['SEED'] = randint(1, 1E9)

    system = cp2k_info['SYSTEM']
    if system in ['PBC_CRYSTAL', 'NEUTRAL_CRYSTAL']:
        #cp2k_info.update({ 'LIST_ACTIVATED' : task_info['LIST_ACTIVATED']})
        from utils import FSSHOSCrystal as Config
    else:
        sys.exit()

    config = Config(cp2k_info, paths, RESTART=restart_info)
    print "GO FOR RUN %d" % cp2k_info['NDIR']
    config.run(cp2k_info['NDIR'])
    if task_info.get('LIGHT', False):
        shorten_log_file(cp2k_info['NDIR'])

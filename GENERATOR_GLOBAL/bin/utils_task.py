# standard modules
import string, re, struct, sys, math, os, time
import numpy as np
import imp
from collections import deque

# custom modules
from utils import *



# FUNCTION TO RUN CP2K
def run_fssh_from_diabat(cp2k_info, task_info, paths):
    # DETERMINE INITIAL STATE TO GET BOLTZMANN RATIO
    scaling = cp2k_info['SCALING']
    initial = Dir('initial/config-%s' % (task_info['FILE_INIT']), paths)



    cp2k_info.update({'STEPS' : int(task_info['LENGTH_FS'] / cp2k_info['TIMESTEP'] ) } )
    cp2k_info.update({'PRINT':  int(task_info['LENGTH_FS'] / cp2k_info['TIMESTEP']) } )
    cp2k_info.update({'PRINT_FSSH': int( 1 / cp2k_info['TIMESTEP']) })

    initial.checkdir()
    paths.update({'initial': initial.path})
    systems = InputFile( initial.path + '/system.info').dict
    cp2k_info.update(systems)

    system = cp2k_info['SYSTEM']
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    config = Config(cp2k_info, paths)
    ndir = cp2k_info['NDIR']
    print "GO FOR RUN %d" % ndir
    config.run(ndir)
    if task_info.get('LIGHT', False):
        with open('run-%s/run.log' % ndir) as oldlog, open('run-%s/new.log' % ndir, 'w' ) as newlog:
            result = ' '.join([next(oldlog) for x in xrange(100)])
            result += ' '.join(deque(oldlog, 100))
            newlog.write(result)
        os.system('mv run-%s/new.log run-%s/run.log' % (ndir, ndir))
        os.system('rm run-%s/run-r-1.out' % ndir)
        os.system('rm run-%s/run-mix-1.ener' % ndir)



def run_equilibrium_fssh(cp2k_info, task_info, paths, dict_equilibrium):
    # DETERMINE INITIAL STATE TO GET BOLTZMANN RATIO
    scaling = cp2k_info['SCALING']
    population = [ task_info.get('NUMBER_CONFIG') * task_info['NUMBER_REPEAT'] ]
    for state in range(1, task_info['NUMBER_ADIABAT']):
        pop = max( round(  dict_equilibrium[scaling][state]  * task_info.get('NUMBER_CONFIG') * task_info['NUMBER_REPEAT']) , 1)
        population.append(pop)
        population[0] = population[0] - pop
    indice = (cp2k_info['REPEAT'] + 1) + (cp2k_info['INIT'] - 1) * task_info['NUMBER_REPEAT']
    this_pop = 0
    for state in range(len(population)):
        this_pop += population[state]
        if indice <= this_pop:
            initial = Dir('initial/state-%s-scaling-%s-%s' % (state + 1, scaling,  task_info['FILE_INIT']), paths)
            cp2k_info.update({'FIRST_ADIABAT': state + 1})
            break
    print population, indice, cp2k_info['FIRST_ADIABAT']

    cp2k_info.update({'STEPS' : int(task_info['LENGTH_FS'] / cp2k_info['TIMESTEP'] ) } )
    cp2k_info.update({'PRINT':  int(task_info['LENGTH_FS'] / cp2k_info['TIMESTEP']) } )
    cp2k_info.update({'PRINT_FSSH': int( 1 / cp2k_info['TIMESTEP']) })

    initial.checkdir()
    paths.update({'initial': initial.path})
    systems = InputFile( initial.path + '/system.info').dict
    cp2k_info.update(systems)

    system = cp2k_info['SYSTEM']
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    config = Config(cp2k_info, paths)
    ndir = cp2k_info['NDIR']
    print "GO FOR RUN %d" % ndir
    config.run(ndir)
    if task_info.get('LIGHT', False):
        with open('run-%s/run.log' % ndir) as oldlog, open('run-%s/new.log' % ndir, 'w' ) as newlog:
            result = ' '.join([next(oldlog) for x in xrange(100)])
            result += ' '.join(deque(oldlog, 100))
            newlog.write(result)
        os.system('mv run-%s/new.log run-%s/run.log' % (ndir, ndir))
        os.system('rm run-%s/run-r-1.out' % ndir)
        os.system('rm run-%s/run-mix-1.ener' % ndir)


def run_sample_bo( cp2k_info, task_info, paths):
    initial_name =  task_info['FILE_INIT'][cp2k_info['FIRST_ADIABAT']]
    initial = Dir('initial/' +  initial_name, paths)
    initial.checkdir()
    paths.update({'initial': initial.path})
    systems = InputFile( initial.path + '/system.info').dict
    cp2k_info.update(systems)

    cp2k_info.update({'STEPS' : int(task_info['LENGTH_FS'] / cp2k_info['TIMESTEP'] ) } )
    cp2k_info.update({'PRINT':  int(cp2k_info['STEPS']/ task_info['OUTPUT_CONFIG'] ) } )
    cp2k_info.update({'PRINT_FSSH': int( 1 / cp2k_info['TIMESTEP']) })

    system = cp2k_info['SYSTEM']
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    output = Dir('output/state-%s-scaling-%s-%s' % (cp2k_info['FIRST_ADIABAT'], cp2k_info['SCALING'], paths['bucket'].split('/')[-1]) , paths )
    output.rm_mkdir()

    config = Config(cp2k_info, paths)
    ndir = cp2k_info['NDIR']
    print "GO FOR RUN %d" % ndir
    config.run(ndir)

    if os.path.exists('run-%d' % (ndir)):
        fssh_parcel = FSSHParcelBO(cp2k_info, paths, output)
        fssh_parcel.gather_vel_coord(config.ndir)
        fssh_parcel.create_system_info(output_path=output.path)


def run_fssh( cp2k_info):
    task_info = cp2k_info.get('INPUTS_DICT')
    paths  = cp2k_info.get('PATHS_DICT')

    if task_info['INITIALIZATION'] == 'ADIABATIC':
        dict_equilibrium = cp2k_info['DICT_EQUILIBRIUM']
        run_equilibrium_fssh(cp2k_info, task_info, paths, dict_equilibrium)
    elif task_info['INITIALIZATION'] == 'SAMPLE_BO':
        run_sample_bo(cp2k_info, task_info, paths)
    elif task_info['INITIALIZATION'] == 'DIABATIC':
        run_fssh_from_diabat(cp2k_info, task_info, paths)
    else:
        print "NO METHOD IMPLEMENTED FOR INITIALIZATION"
        sys.exit()




#def estimate_boltzmann_ratio(scaling, paths, nconfig):  # scaling in Ha
#    scalings = np.array(jacob_coupling.keys()) * 0.0000367493  # convert in Ha
#    couplings = np.array(jacob_coupling.values()) * 0.001  # convert in eV
#    fit = np.polyfit(scalings, couplings, 1)
#    coupling =  fit[0] * scaling + fit[1]
                 #
#    marcus = imp.load_source('scripts', paths.get('bucket') + '/scripts/marcus.py')
#    temperature = 300
#    reorga = 0.100
#    free_energy = 0
                 #
#    return np.rint( nconfig * marcus.calculate_boltzman_ratio(reorga, free_energy, coupling, temperature, state = 'Ground') )
                 #
                 #
                 #
#def estimate_scaling(coupling):  # coupling in eV
#    #scalings = np.array(jacob_coupling.keys()) * 0.0000367493  # convert in Ha
#    scalings = np.array(cp2k_coupling.keys()) # scaling in Ha
#    couplings = np.array(cp2k_coupling.values()) * 27.211399 # convert in eV
#    fit = np.polyfit(scalings, couplings, 1)
#    scaling = (coupling - fit[1]) / fit[0]
#    #coupling =  fit[0] * scaling + fit[1]
#    return scaling # in Ha
                 #
                 #
                 #
                 #




# FUNCTION TO FIND CORRESPONDING ENERGIES
def check_temperature_diff(temp1, temp2, target):
    criteria = 3  # Kelvin
    if ((abs(temp1 - target)) < criteria) and (abs(abs(temp1 - temp2) - 5) < criteria) and (temp2 < temp1):
        # print temp1, temp2
        return True
    else:
        # print temp1, temp2
        return False


def compare_last_energy(run_results, target):
    with open('run-%s/run-1.ener' % run_results[0]) as file1, open('run-%s/run-1.ener' % run_results[1]) as file2:
        temp1 = float(list(deque(file1, 1))[0].split()[3])
        temp2 = float(list(deque(file2, 1))[0].split()[3])
        return check_temperature_diff(temp1, temp2, target)


def compare_all_energies(run_results, target):
    lines = 1
    with open('run-%s/run-1.ener' % run_results[0]) as file1, open('run-%s/run-1.ener' % run_results[1]) as file2:
        for line1, line2 in zip(file1.readlines()[1:], file2.readlines()[1:]):
            temp1 = float(line1.split()[3])
            temp2 = float(line2.split()[3])
            if check_temperature_diff(temp1, temp2, target):
                return lines
            else:
                lines += 1
    return -1


def clean_ante_run(run_results):
    for ndir in run_results:
        antedir = ndir - 2
        if (os.path.exists('run-%d' % antedir)) and (antedir != 0) and (antedir != 1):
            os.system('rm run-%d/run-*.xyz' % antedir)




# FUNCTION TO RUN FIST
def run_fist_long_nvt(cp2k_info):
    inputs = cp2k_info.get('INPUTS_DICT')
    inputs.update(cp2k_info)
    paths  = cp2k_info.get('PATHS_DICT')
    system = inputs.get('SYSTEM')
    cp2k_infostate_temp = cp2k_info.get('DICT_STATE_TEMP')
    inputs.update(cp2k_infostate_temp[cp2k_info['FIRST_ADIABAT']])

    print "1. CONSTRUCT THE ORGANIC CRYSTAL."
    if system == 'CRYSTAL':
        from utils import OSCluster as Structure
        from utils import CP2KOS as Config
    elif system == 'SOLVENT':
        from utils import OSwSolvent as Structure
        from utils import CP2KOSwSolvent as Config
    else:
        sys.exit()

    ndir = cp2k_info['NDIR']
    output = Dir('output/state-%s-density-%s-%s-run-%s' % (cp2k_info['FIRST_ADIABAT'], cp2k_info['DENSITY'], paths['bucket'].split('/')[-1], ndir ) , paths)
    output.rm_mkdir()
    paths.update({ 'output' : output.path })

    structure = Structure(inputs, paths, ndir)

    print "2. RUN CP2K"
    print "GO FOR RUN %d" % ndir
    config_nvt = Config(inputs, paths, ENSEMBLE='NVT', STEPS=inputs['NEQ'], RESTART=None, TEMPLATE_FILE='FIST_ETHYLENE_CONSTRAINT.template')
    ndir = config_nvt.run(ndir)

    return config_nvt.ndir


def run_fist_short_nvt(cp2k_info):
    inputs = cp2k_info.get('INPUTS_DICT')
    inputs.update(cp2k_info)
    paths  = cp2k_info.get('PATHS_DICT')
    system = inputs.get('SYSTEM')
    cp2k_infostate_temp = cp2k_info.get('DICT_STATE_TEMP')
    inputs.update(cp2k_infostate_temp[cp2k_info['FIRST_ADIABAT']])

    inputs['NCONFIG'] = inputs['NPROD']

    if system == 'CRYSTAL':
        from utils import CP2KOS as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolvent as Config
    else:
        sys.exit()

    ndir = cp2k_info['NDIR']
    print "2. RUN CP2K"
    print "GO FOR RUN %d" % ndir
    config_nvt = Config(inputs, paths, ENSEMBLE='NVT', STEPS=inputs['NPROD'], RESTART=cp2k_info['RESTART'], TEMPLATE_FILE='FIST_ETHYLENE_CONSTRAINT_with_vel')
    ndir = config_nvt.run(ndir)

    return config_nvt.ndir

def run_fist_nvt_nve_extract(dict_):
    inputs = dict_.get('INPUTS_DICT')
    inputs.update(dict_)
    paths  = dict_.get('PATHS_DICT')
    system = inputs.get('SYSTEM')
    #dict_state_temp = dict_.get('DICT_STATE_TEMP')
    #inputs.update(dict_state_temp[dict_['FIRST_ADIABAT']])

    print "1. CONSTRUCT THE ORGANIC CRYSTAL."
    if system == 'CRYSTAL':
        from utils import OSCluster as Structure
        from utils import CP2KOS as Config
    elif system == 'SOLVENT':
        from utils import OSwSolvent as Structure
        from utils import CP2KOSwSolvent as Config
    else:
        sys.exit()

    if dict_.get('DENSITY'):
        output = Dir('output/config-%s-%s' % ( dict_['DENSITY'], paths['bucket'].split('/')[-1]) , paths )
    else:
        output = Dir('output/config-%s' % ( paths['bucket'].split('/')[-1]) , paths )
    output.rm_mkdir()
    paths.update({ 'output' : output.path })

    structure = Structure(inputs, paths)

    ndir = dict_['NDIR']
    print "2. RUN CP2K"
    print "GO FOR RUN %d" % ndir
    config_nvt = Config(inputs, paths, ENSEMBLE='NVT', STEPS=inputs['NEQ'], RESTART=None, TEMPLATE_FILE='FIST_without_vel.template')
    ndir = config_nvt.run(ndir)

    config_nve = Config(inputs, paths, ENSEMBLE='NVE', STEPS=inputs['NPROD'], RESTART=config_nvt.ndir, TEMPLATE_FILE='FIST_TEMPLATE')
    ndir = config_nve.run(ndir)

    if os.path.exists('run-%d' % (config_nve.ndir)):
        fssh_parcel = FSSHParcel(inputs, paths)
        fssh_parcel.gather_vel_coord(config_nve.ndir, output_path=output.path)
        fssh_parcel.create_system_info(output_path=output.path)




# FUNCTION TO EXTRACT CONFIG
def extract_all_config(run_results, inputs, paths):
    for ndir in run_results:
        if os.path.exists('run-%d' % (ndir)):
            output = Dir('output/state-%s-density-%s-%s-run-%s' % (inputs['FIRST_ADIABAT'], inputs['DENSITY'], paths['bucket'].split('/')[-1], ndir), paths)
            output.rm_mkdir()
            fssh_parcel = FSSHParcel(inputs, paths)
            fssh_parcel.gather_vel_coord(ndir, output_path=output.path)
            fssh_parcel.create_system_info(output_path=output.path)


def extract_one_config(run_results, inputs, paths, line):
    for ndir in run_results:
        if os.path.exists('run-%d' % (ndir)):
            output = Dir('output/run-%s-density-%s-%s' % (ndir, inputs['DENSITY'], paths['bucket'].split('/')[-1]), paths)
            output.rm_mkdir()
            fssh_parcel = FSSHParcel(inputs, paths)
            fssh_parcel.gather_vel_coord(ndir, output_path=output.path, line=line)
            fssh_parcel.create_system_info(output_path=output.path)




## DATE FROM PREVIOUS RUNS
#jacob_coupling = {  # meV : meV
#    89: 8,
#    133: 11,
#    177: 15,
#    266: 23,
#    355: 30,
#    532: 45,
#    710: 61,
#    887: 76,
#    1330: 114,
#    1774: 151,
#    2217: 189,
#    2661: 227,
#    3104: 265,
#    3548: 303,
#    3991: 341,
#    4435: 378,
#    4878: 418,
#    5322: 455
#}
 #
#cp2k_coupling = { # Ha : Ha
#    9E-05 : 1.04E-05,
#    0.0005 : 6.05E-05,
#    0.004  : 0.000505,
#    0.02   : 0.0027
#}
 #
## This data are taken from: extract-scaling-adiabat-TASK271-SAMPLE-BO-GROUND-100ps-170516-1dca75020bb2d36bb5841283e6867e5a-1705161812
#dict_1_Simulation = {  \
# 5e-05 : 0.0252513034296,   \
# 0.0001 : 0.0253138065831,   \
# 0.0005 : 0.0253136759442,   \
# 0.001 : 0.0251646663009,   \
# 0.003 : 0.0240682953401,   \
# 0.005 : 0.0223558285774,   \
# 0.008 : 0.0181561059991,   \
# 0.01 : 0.0151159085859,   \
# 0.03 : 0.000904812300885,   \
# 0.05 : 2.80776025514e-05
# }
#
## This data are taken from: extract-scaling-adiabat-TASK271-SAMPLE-BO-30ps-170517-b7905f4556ed79fc2323eeec5105416e-1705172113
#dict_scaling_excited_pop = {  \
# 5e-05 : 52.2374550503,   \
# 0.0001 : 52.4626889155,   \
# 0.0005 : 54.2785150554,   \
# 0.001 : 56.5831263004,   \
# 0.003 : 66.1847344636,   \
# 0.005 : 76.3807141674,   \
# 0.008 : 92.7009484113,   \
# 0.01 : 104.188794639,   \
# 0.02 : 166.8687697,   \
# 0.03 : 234.597679201,   \
# 0.05 : 375.850102623
#}


 #!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy as np
from multiprocessing import Pool, cpu_count

# custom modules
from utils import *


def run_fssh( dict_):
    inputs = dict_.get('INPUTS_DICT')
    paths  = dict_.get('PATHS_DICT')

    initial_name =  'state-%s-density-%s-%s' % (dict_['FIRST_ADIABAT'], dict_['DENSITY'], inputs['FILE_INIT'])
    initial = Dir('initial/' +  initial_name, paths)
    initial.checkdir()
    paths.update({'initial': initial.path})
    systems = InputFile( initial.path + '/system.info').dict
    inputs.update(systems)


    system = inputs['SYSTEM']
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    output = Dir('output/state-%s-scaling-%s-%s' % (dict_['FIRST_ADIABAT'], dict_['SCALING'], paths['bucket'].split('/')[-1]) , paths )
    output.rm_mkdir()

    inputs.update(dict_)
    config = Config(inputs, paths, **dict_)
    ndir = dict_['NDIR']
    print "GO FOR RUN %d" % ndir
    config.run(ndir)

    if os.path.exists('run-%d' % (ndir)):
        fssh_parcel = FSSHParcelBO(inputs, paths, output)
        fssh_parcel.gather_vel_coord(config.ndir)
        fssh_parcel.create_system_info(output_path=output.path)

def estimate_scaling(coupling):  # coupling in eV
    jacob_coupling = {
        89: 8,
        133: 11,
        177: 15,
        266: 23,
        355: 30,
        532: 45,
        710: 61,
        887: 76,
        1330: 114,
        1774: 151,
        2217: 189,
        2661: 227,
        3104: 265,
        3548: 303,
        3991: 341,
        4435: 378,
        4878: 418,
        5322: 455
    }

    scalings = np.array(jacob_coupling.keys()) * 0.0000367493  # convert in Ha
    couplings = np.array(jacob_coupling.values()) * 0.001  # convert in eV
    fit = np.polyfit(scalings, couplings, 1)
    scaling = (coupling - fit[1]) / fit[0]
    #coupling =  fit[0] * scaling + fit[1]

    return scaling

def round_to_1(x):
    return round(x, -int(np.floor(np.log10(np.abs(x)))))

def get_list_scaling(number, reorga):
    division = { 1 : 1, 2 : 2, 3 : 3, 4 :5, 5 : 10, 6 : 100, 7 : 1000}
    scalings = []
    for i in range(1, number + 1):
        div = division[i]
        coupling = reorga / div
        scaling = estimate_scaling(coupling)
        scalings.append(round_to_1(scaling))
    return scalings


def main(inputs, paths):

    
    task = {
        'KIND_RUN'  : 'TASK234-SAMPLE-BO',
        'TEMPLATE_FILE' : 'FSSH_CORE.template',
        'FORCEFIELD_FILE' : 'FSSH_FF.template',
        'NUMBER_INIT'     : 1,
        'NUMBER_RANDOM'   : 1,
        'FILE_INIT' : 'TASK234-SAMPLE-BO-FIRST-ETHYLENE-170404-a5e55a180874b5807f102df3d41df810',
        'STEPS'        : 100,
        'PRINT'        : 1,
        'NUMBER_SCALING' : 7,
        'REORGANIZATION_ENERGY' : 0.1  # eV
            }
    inputs.update(task)

    list_propagation = ['BORN_OPPENHEIMER']
    list_init     = range(1, inputs.get('NUMBER_INIT') + 1)
    list_repeat   = range(inputs.get('NUMBER_RANDOM'))
    list_scaling = get_list_scaling( inputs['NUMBER_SCALING'], inputs['REORGANIZATION_ENERGY']  )
    list_first_adiabat = [1, 2]
    list_density = [0.001]
    list_cc_charged = [1.369]

    mega_list = [ { 'PROPAGATION' : prop,
                    'INIT'              : init ,
                    'REPEAT'             : repeat,
                    'SCALING'            : scaling,
                    'FIRST_ADIABAT'            : adiabat,
                    'DENSITY'             : density,
                    'CC_CHARGED'          : cc_charged
                   }
                  for prop in list_propagation
                  for init      in list_init
                  for repeat    in list_repeat
                  for scaling in list_scaling
                  for adiabat in list_first_adiabat
                  for density in list_density
                  for cc_charged in list_cc_charged
                  ]


    # SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
    bucket = Bucket(inputs)
    bucket.name()
    paths.update({'bucket': bucket.path})

    task = Dir(inputs.get('INPUT_INFO'))
    paths.update( {'task' : task.path} )

    templates = Dir('templates', paths)
    templates.checkdir()
    templates.clean()

    bin = Dir('bin', paths)
    bin.checkdir()

    output_total = Dir('output')
    output_total.rm_mkdir()


    # PREPARE THE MEGA_LIST FOR POOL
    for ndir in range(len(mega_list)):
        mega_list[ndir].update({ 'NDIR' : ndir,
                                 'PATHS_DICT' : paths,
                                 'INPUTS_DICT' : inputs
                                 })


    # RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
    nworker = inputs['NWORKER']
    if nworker == -1:
        nworker = cpu_count()

    if nworker == 0:
        for dict_ in mega_list:
            run_fssh(dict_)
    else:
        pool = Pool(nworker)
        pool.map( run_fssh, mega_list)




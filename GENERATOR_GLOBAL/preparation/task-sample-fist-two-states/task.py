#!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy
from multiprocessing import Pool, cpu_count

# custom modules
from utils import *
from utils_task import *





def main(inputs, paths):
    print """
    """

    task = {
        'KIND_RUN'   : 'TASK271-SAMPLE-FIST',
        'NCONFIG' : 1,
        'NEQ'  : 10000,
        'NPROD'  : 10000,
        'SYSTEM': 'SOLVENT',
        'MOL_NAME'      : 'ETHYLENE',
        'FILEMOL'       : 'ethylene.coord',
        'FILECRYSTAL'   : 'ethylene_dimer.coord',
        'NAME_SOLVENT': 'NE',
        'SOLVENT': 'Ne',
        'PERIODIC' : 'XYZ',
        'SIZE_BOX'     : [60.0, 60.0, 60.0],
        'CLOSEST_DIST'  : 5,
        'VECTA'        : [3.527, 0.784, -0.166],
        'VECTB'        : [0.0, 0.0, 0.0],
        'VECTC'        : [0.0, 0.0, 0.0],
        'SIZE_CRYSTAL' : [2, 1, 1],
        'COORD_CHARGE' : [2, 1, 1],
        'RCUT'         : 12,
        'TEMPERATURE'   : 300
    }
    inputs.update(task)

    list_density     = [0.001]
    list_state        = [1, 2]
    dict_state_temp   = {
        1  : {
            'CC_CHARGED'  : 1.369
        },
        2  : {
            'CC_CHARGED'  : 1.3239,
        }
    }
    mega_list = [ {
                    'DENSITY'          : density,
                    'FIRST_ADIABAT'    : state
                   }
                  for density  in list_density
                  for state    in list_state
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
                                 'INPUTS_DICT' : inputs,
                                 'DICT_STATE_TEMP' : dict_state_temp
                                 })


    # RUN THE FIRST NVT CALCULATION
    # IN SERIE
    nworker = 0
    if nworker == 0:
        run_results = []
        for cp2k_info in mega_list:
            run_results.append(run_fist_long_nvt(cp2k_info))
    else:
        from multiprocessing import Pool, cpu_count
        if nworker == -1:
            nworker = cpu_count()
        pool = Pool(nworker)
        run_results = pool.map( run_fist_long_nvt, mega_list)


    # CHECK IF LAST CONFIG GET THE CRITERIA
    # IF NOT, REDO SMALLER NVT RUNS UNTIL ONE CONFIG GET THE CRITERIA
    target = 300
    converged_ = compare_last_energy(run_results, target)
    line_config = inputs['NEQ']
    while not converged_:
        for ndir in range(len(mega_list)):
            mega_list[ndir].update({'NDIR': run_results[ndir] + 2,
                                    'RESTART' : run_results[ndir]
                                    })
        run_results = []
        for cp2k_info in mega_list:
            run_results.append(run_fist_short_nvt(cp2k_info))
        line_converged = compare_all_energies(run_results, target)
        if line_converged == -1:
            converged_ = False
            clean_ante_run(run_results)
        else:
            converged_ = True
            line_config = line_converged


    # EXTRACT THE CONFIG WHICG GET THE CRITERIA
    extract_one_config(run_results, inputs, paths, line_config)


    # COPY ANALYSER
    os.system('cp %sanalyser.py %s' % (paths.get('task'), paths.get('bucket')))

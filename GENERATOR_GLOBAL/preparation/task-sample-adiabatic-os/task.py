
 #!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy
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
            }
    inputs.update(task)

    list_propagation = ['BORN_OPPENHEIMER']
    list_init     = range(1, inputs.get('NUMBER_INIT') + 1)
    list_repeat   = range(inputs.get('NUMBER_RANDOM'))
    list_scaling = [0.0001, 0.001, 0.01, 0.03, 0.05, 0.07, 0.1]
    #list_scaling = [0.0001]
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




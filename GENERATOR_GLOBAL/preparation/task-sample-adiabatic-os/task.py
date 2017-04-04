
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
    system = inputs.get('SYSTEM')
    state_dict = dict_.get('STATE_DICT')
    if system == 'CRYSTAL':
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    initial = Dir('initial/' + dict_['INITIAL'], paths)
    initial.checkdir()
    paths.update({'initial': initial.path})

    output = Dir('output/%s-scaling-%s-%s' % (dict_['INITIAL'], dict_['SCALING'], paths['bucket'].split('/')[-1]) , paths )
    output.rm_mkdir()

    dict_.update(
        {'FIRST_DIABAT' : state_dict[ dict_['INITIAL'] ]}
    )

    config = Config(inputs, paths, **dict_)
    ndir = dict_['NDIR']
    print "GO FOR RUN %d" % ndir
    config.run(ndir)

    if os.path.exists('run-%d' % (ndir)):
        fssh_parcel = FSSHParcelBO(inputs, paths, output)
        fssh_parcel.gather_vel_coord(config.ndir)




def main(inputs, paths):

    
    task = {
        'KIND_RUN'  : 'TASK234-SAMPLE-BO',
        'TEMPLATE_FILE' : 'FSSH_CORE.template',
        'FORCEFIELD_FILE' : 'FSSH_FF.template',
        'FILE_INIT' : 'initial',
        'PERIODIC' : 'XYZ',
        'NUMBER_INIT' : 1,
        'NUMBER_RANDOM' : 1,
        'SYSTEM' : 'SOLVENT',
        'MOL_NAME' :     'ETHYLENE',
        'NAME_SOLVENT' : 'NE',
        'SOLVENT'      : 'Ne',
        'NATOMS'       : 136,
        'NATOM_MOL'    : 6,
        'SIZE_BOX'     : [60.0, 60.0, 60.0],
        'VECTA'        : [3.527, 0.784, -0.166],
        'VECTB'        : [0.0, 0.0, 0.0],
        'VECTC'        : [0.0, 0.0, 0.0],
        'SIZE_CRYSTAL' : [2, 1, 1],
        'COORD_CHARGE' : [2, 1, 1],
        'STEPS'        : 100,
        'PRINT'        : 1,
        'RCUT'      :    12,
        'CC_CHARGED'   : 1.369
            }
    inputs.update(task)

    list_propagation = ['BORN_OPPENHEIMER']
    list_init     = range(1, inputs.get('NUMBER_INIT') + 1)
    list_repeat   = range(inputs.get('NUMBER_RANDOM'))
    list_scaling = [0.0001, 0.001, 0.01, 0.03, 0.05, 0.07, 0.1]
    #list_scaling = [0.0001]
    state_dict = {
        'ground-state': 1,
        'excited-state': 2
    }


    mega_list = [ { 'PROPAGATION' : prop,
                    'INIT'              : init ,
                    'REPEAT'             : repeat,
                    'SCALING'            : scaling,
                    'INITIAL'            : initial
                   }
                  for prop in list_propagation
                  for init      in list_init
                  for repeat    in list_repeat
                  for scaling in list_scaling
                  for initial in state_dict.keys()
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

    for element in state_dict.keys():
        os.system(' cp -r %s/%s %s' % (paths.get('task'), element, paths.get('bucket')))


    # PREPARE THE MEGA_LIST FOR POOL
    for ndir in range(len(mega_list)):
        mega_list[ndir].update({ 'NDIR' : ndir,
                                 'PATHS_DICT' : paths,
                                 'INPUTS_DICT' : inputs,
                                 'STATE_DICT' : state_dict
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




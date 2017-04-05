#!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy

# custom modules
from utils import Dir, FSSHParcel, Bucket


def run_fist(dict_):
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


def main(inputs, paths):
    print """
    START THE DIABAT OF SAMPLING OF AN ORGANIC CRYSTAL 
    """

    # SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
    bucket = Bucket(inputs)
    bucket.name()
    paths.update({'bucket': bucket.path})

    templates = Dir('templates', paths)
    templates.checkdir()
    templates.clean()

    bin = Dir('bin', paths)
    bin.checkdir()

    output = Dir('output', paths)
    output.rm_mkdir()


    # PREPARE THE MEGA_LIST FOR POOL
    mega_list = [{}]
    for ndir in range(len(mega_list)):
        mega_list[ndir].update({ 'NDIR' : ndir*2,
                                 'PATHS_DICT' : paths,
                                 'INPUTS_DICT' : inputs
                                 })


    # RUN THE CALCULATIONS, SERIE OR PARALLEL ACCORDING TO THE NWORKER VARIABLE
    nworker = 0
    #nworker = inputs['NWORKER']
    #if nworker == -1:
    #    nworker = cpu_count()

    if nworker == 0:
        for dict_ in mega_list:
            run_fist(dict_)
    #else:
    #    pool = Pool(nworker)
    #     pool.map(run_fist, mega_list)


    # COPY ANALYSER
    os.system('cp %sanalyser.py %s' % (paths.get('task'), paths.get('bucket')))
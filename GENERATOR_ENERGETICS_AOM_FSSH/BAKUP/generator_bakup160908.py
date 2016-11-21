#!/usr/bin/python

import string, re, struct, sys, math, os, time

def generate_coord_files(FILEIN, FILEOUT, atom_number, atom_mol, diabat_number):
    template = open(FILEIN,'r')
    lines    = template.readlines()
    template.close
    atom_label = []
    atom_coord = []
    for atom_index in range(atom_number):
        l = string.strip(lines[atom_index])
        info = re.split('\s+',l)
        atom_label.append(info[0])
        atom_coord.append([info[1], info[2], info[3]])
    for diabat_index in range(diabat_number):
        fileout = open( FILEOUT + '_' + str(diabat_index+1) + '.init', 'w')
        for atom_index in range(atom_number):
            index = ( atom_index + diabat_index * atom_mol) % (atom_number)
            result = '%s  %s %s %s\n' \
                     %( atom_label[index], \
                        atom_coord[atom_index][0], \
                        atom_coord[atom_index][1], \
                        atom_coord[atom_index][2] )
            fileout.write(result)
    fileout.close()
 
def generate_input_file_dimer(METHOD, SURFACE = 1, SCALING = 0.0065, TIMESTEP = 0.5, STEPS = 10, PRINT = 5 ):
   cmd = 'sed                     \
        -e \'s/sedMETHOD/%s/g\'   \
        -e \'s/sedSURFACE/%s/g\'  \
        -e \'s/sedSCALING/%s/g\'  \
        -e \'s/sedTIMESTEP/%s/g\' \
        -e \'s/sedSTEPS/%s/g\'    \
        -e \'s/sedPRINT/%s/g\'    \
       DIMER_ETHYLENE_TEMPLATE > run.inp' \
      % (METHOD, SURFACE, SCALING, TIMESTEP, STEPS, PRINT) 
   os.system(cmd)
   generate_coord_files(FILEIN='DIMER_COORD',FILEOUT='COORD',atom_number=12, atom_mol=6,diabat_number=2)



try : 
         DATE         = sys.argv[1]
         CP2K_VERSION = sys.argv[2]
except :
         raise SystemExit

MODULEPATH = '/scratch/grudorff/modulefiles module load mpich'
PWD = os.getcwd()
PATH_EXE = '/scratch/grudorff/antoine/bin'
PATH_TEMPLATE = PWD  + '/TEMPLATE/'
PATH_COORD = PWD + '/COORD/'
exe_cp2k = PATH_EXE + '/cp2k_nonadiabatic_scutum_160427.sopt'

BUCKET='ENERGIES_TEST_%s_%s' % (DATE, CP2K_VERSION)
os.chdir('..')
test_bucket = os.path.exists(BUCKET)
if (test_bucket):
   print "THE DIRECTORY EXISTS ALREADY !"
   sys.exit()
else:
   os.mkdir(BUCKET)
   os.chdir(BUCKET)

number_dir = -1


print """!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
STUDY BORN_OPPENHEIMER FOR DIFFERENT SCALING FACTOR (COUPLING)
"""
scalingList = ['0.0065', '0.0130', '0.0325', '0.0652', '0.1956', '0.3912', '0.9780']  # in eV
stateList   = ['1', '2']

CONFIG = 1
METHOD = 'BORN_OPPENHEIMER'
for state in range(2):
    statep = state + 1

    for scalingInd in range(len(scalingList)):
        scaling = scalingList[scalingInd] 
        number_dir = number_dir + 1
        run = 'run%d' % number_dir
        number_config = 1
        name_sys = 'BO%d_C%f'
        print """
        RUN = %s
        METHOD  = %s
        SURFACE = %d
        CONFIG  = %d
        SCALING FACTOR = %f
        """ % (run, METHOD, state, CONFIG, float(scaling))
        os.mkdir(run)
        os.chdir(run)
        cmd = 'cp %s/* .' % PATH_TEMPLATE
        os.system(cmd)
        cmd = 'cp %s/COORD_%d COORD.init' % ( PATH_COORD, CONFIG) 
        generate_input_file_dimer(METHOD=METHOD, SURFACE=state, SCALING=scaling)
        cmd = exe_cp2k + ' run.inp > run.out 2>   run.err'
        os.system(cmd)
        os.chdir('..')





print """
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

STUDY SURFACE HOPPING FOR DIFFERENT SCALING FACTOR

"""
CONFIG = 1
METHOD = 'FSSH'
scalingList = ['0.0065', '0.0130', '0.0325', '0.0652', '0.1956', '0.3912', '0.9780']  # in eV
for scalingInd in range(len(scalingList)):
    scaling = scalingList[scalingInd]
    number_dir = number_dir + 1
    run = 'run%d' % number_dir
    number_config = 1
    name_sys = 'BO%d_C%f'
    print """
    RUN = %s
    METHOD  = %s
    CONFIG  = %d
    SCALING FACTOR = %f
    """ % (run, METHOD, CONFIG, float(scaling))
    os.mkdir(run)
    os.chdir(run)
#        subdir = "CONFIG_" + config
#        os.mkdir(subdir)
#        os.chdir(subdir)
    cmd = 'cp %s/* .' % PATH_TEMPLATE
    os.system(cmd)
    cmd = 'cp %s/COORD_%d COORD.init' % ( PATH_COORD, CONFIG)
#            generate_input_file()
#            cmd = '%s %s.inp > %s.out 2>   %s.out'
#            os.system(cmd)
    os.chdir('..')






print """
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

STUDY SURFACE HOPPING FOR DIFFERENT TIMESTEP 

"""
timestepList = ['0.05', '0.1', '0.5', '1']
scaling = 0.0065
for timestepInd in range(len(timestepList)):
    timestep = timestepList[timestepInd]
    number_dir = number_dir + 1
    run = 'run%d' % number_dir
    number_config = 1
    name_sys = 'BO%d_C%f'
    print """
    RUN = %s
    METHOD  = %s
    CONFIG  = %d
    TIMESTEP = %f
    """ % (run, METHOD, CONFIG, float(timestep))
    os.mkdir(run)
    os.chdir(run)
#        subdir = "CONFIG_" + config
#        os.mkdir(subdir)
#        os.chdir(subdir)
    cmd = 'cp %s/* .' % PATH_TEMPLATE
    os.system(cmd)
    cmd = 'cp %s/COORD_%d COORD.init' % ( PATH_COORD, CONFIG)
#            generate_input_file()
#            cmd = '%s %s.inp > %s.out 2>   %s.out'
#            os.system(cmd)
    os.chdir('..')    









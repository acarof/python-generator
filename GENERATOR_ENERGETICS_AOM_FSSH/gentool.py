
import string, re, struct, sys, math, os, time
import hashlib
import subprocess


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


 
def generate_input_file_dimer(METHOD, FILECOORD, FILEVELOC,DIABAT = 1, SCALING = 0.0065, TIMESTEP = 0.5, STEPS = 10, PRINT = 5 ):
   cmd = 'sed                     \
        -e \'s/sedMETHOD/%s/g\'   \
        -e \'s/sedDIABAT/%s/g\'  \
        -e \'s/sedSCALING/%s/g\'  \
        -e \'s/sedTIMESTEP/%s/g\' \
        -e \'s/sedSTEPS/%s/g\'    \
        -e \'s/sedPRINT/%s/g\'    \
       DIMER_ETHYLENE_TEMPLATE > run.inp' \
      % (METHOD, DIABAT, SCALING, TIMESTEP, STEPS, PRINT) 
   os.system(cmd)
   generate_coord_files(FILEIN=FILECOORD, FILEOUT='COORD', atom_number=12, atom_mol=6,diabat_number=2)




def create_and_run(CONFIGLIST = [1], METHODLIST = ['FSSH'], DIABATLIST = [1], SCALINGLIST  = [0.0065],  \
                   TIMESTEPLIST = [0.5], number_dir = -1, STEPSLIST = [1], NAME_COORD = 'COORD',              \
                   NAME_VELOC = 'VELOC.init'):
    for config in CONFIGLIST:
        for method in METHODLIST:
            for diabat in DIABATLIST:
                for scaling in SCALINGLIST:
                    for timestepInd in range(len(TIMESTEPLIST)):
                        timestep = TIMESTEPLIST[timestepInd]
                        steps = STEPSLIST[timestepInd]
                        number_dir = number_dir + 1
                        run = 'run%d' % number_dir
                        print """
                          RUN = %s
                          CONFIG  = %d
                          METHOD  = %s
                          DIABAT = %d
                          SCALING FACTOR = %f
                          TIMESTEP, LENGTH = %f, %d""" \
                          % (run, config, method, diabat, scaling, timestep, steps)
                        os.mkdir(run)
                        os.chdir(run)
                        os.system('cp %s/* .' % PATH_TEMPLATE)
                        name_coord_config = NAME_COORD + '_' + str(config)
                        name_veloc_config = NAME_VELOC  # TO MODIFY IN THE FUTURE
                        os.system('cp %s/%s    . ' % ( PATH_CONFIG, name_coord_config))
                        os.system('cp %s/%s    . ' % ( PATH_CONFIG, name_veloc_config))  
                        generate_input_file_dimer(METHOD=method, DIABAT=diabat + 1, \
                              SCALING=scaling, TIMESTEP = timestep, STEPS = steps,  \
                              FILECOORD = name_coord_config, FILEVELOC = name_veloc_config)
                        os.system(exe_cp2k + ' run.inp > run.out 2>   run.err')
                        os.chdir('..')
    return number_dir; 



def create_bucket(BUCKET):
    complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
    md5 = hashlib.md5()
    md5.update(complete_time)
    md5name = md5.hexdigest()

    short_time = time.strftime("%y%m%d", time.localtime())

    proc = subprocess.Popen([exe_cp2k ,'--version'], stdout=subprocess.PIPE)
    cp2k_version = proc.communicate()[0]

    for item in cp2k_version.split("\n"):
        if "git:" in item:
            gitline = item.strip()
 
    search_version = re.search('git:(.+?)\n', cp2k_version )
    if search_version:
       cp2k_version = search_version.group(1)

    bucket_name = BUCKET + '_' + short_time + '_' + cp2k_version + '_' + md5name
    os.chdir('..')
    test_bucket = os.path.exists(bucket_name)
    if (test_bucket):
       print "THE DIRECTORY EXISTS ALREADY!"
       sys.exit()
    else:
       os.mkdir(bucket_name)
       os.chdir(bucket_name)


MODULEPATH = '/scratch/grudorff/modulefiles module load mpich'
PWD = os.getcwd()
PATH_EXE = '/scratch/grudorff/antoine/bin'
PATH_TEMPLATE = PWD  + '/TEMPLATE/'
PATH_CONFIG = PWD + '/CONFIG_INITIAL'
exe_cp2k = PATH_EXE + '/cp2k_nonadiabatic_scutum_160427.sopt'












import string, re, struct, sys, math, os, time
import hashlib
import subprocess
from numpy import dot, array


def choose_atom_label(atom, i3d = [-1, -1, -1], i3dcharge = [-1, -1, -1], imol = -1, icharge = -1):
    if (imol == -1) and (icharge == -1):
       if (atom == 'C'):
          if (i3d == i3dcharge):
             atom_label = 'CP'
          else:
             atom_label = 'CN'
       else:
          atom_label = atom
    elif (imol != -1) and (icharge != -1):
       if (atom == 'C'):
          if (imol == icharge):
             atom_label = 'CP'
          else:
             atom_label = 'CN'
       else:
          atom_label = atom
    else:
       print "THERE IS A PROBLEM in choose_atom_label"
    return atom_label


def construct_organic_crystal(MOLFILENAME, MESHA, MESHB, MESHC, SIZE3D, INITCHARGE, CRYSTALFILENAME):
    molfile = open(MOLFILENAME, 'r')
    crystalfile = open(CRYSTALFILENAME, 'w')
    lines   = molfile.readlines()
    molfile.close
    for mol0_index in range(SIZE3D[0]):
        for mol1_index in range(SIZE3D[1]):
            for mol2_index in range(SIZE3D[2]):
                index3d = [mol0_index+1, mol1_index+1, mol2_index+1] 
                for line_index in range(len(lines)):
                    l = string.strip(lines[line_index])
                    info = re.split('\s+',l)
                    atom_label = choose_atom_label(info[0], i3d = index3d, i3dcharge = INITCHARGE)
                    atom_coord = [float(info[1]), float(info[2]), float(info[3])]
                    vec_shift = mol0_index * array(MESHA) + mol1_index * array(MESHB) + mol2_index * array(MESHC)
                    result = '%s  %s\n' \
                    % (atom_label, str(atom_coord + vec_shift).strip('[]') )
                    crystalfile.write(result)
    crystalfile.close


def gather_config(molecule, nmol, natom_mol, icoord_charge, size_crystal, nsteps_eq, nsteps_prod, print_frq):
    icharge    = (float(icoord_charge[0]) - 1 )*float(size_crystal[1])*float(size_crystal[2]) + \
                 (float(icoord_charge[1]) - 1 )*float(size_crystal[2]) + \
                  float(icoord_charge[2]) 
    for iprop in ['vel','pos']:
        #filein = open(molecule + '_CRYSTAL-' + iprop +'-1.xyz', 'r')
        filein = open('run-' + iprop +'-1.xyz', 'r')
        lines = filein.readlines()
        iconfig = 0
        for istep in range(nsteps_eq, nsteps_eq + nsteps_prod, print_frq):
            iconfig = iconfig + 1
            fileout = open( iprop + '-' + str(iconfig), 'w')
            for imol in range(nmol):
                   for iatom in range(natom_mol):
                       index = (istep/print_frq)*(nmol*natom_mol+2) \
                             + 2                    \
                             + imol*natom_mol           \
                             + iatom
                       l = string.strip(lines[index])
                       info = re.split('\s+',l)
                       atom_label = choose_atom_label(info[0], imol = imol, icharge = icharge)
                       atom_xyz = [float(info[1]), float(info[2]), float(info[3])]
                       result = '%s  %s\n' \
                           % (atom_label, str(atom_xyz).strip('[]') )
                       fileout.write(result)
            fileout.close()
                 


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


def input_fist(molecule, timestep, nsteps, print_frq, natoms, sizebox):
    rcut = max(sizebox)/2
#-e \'s/sedPROJECT_NAME/%s_CRYSTAL/g\'  \
    cmd = 'sed \
        -e \'s/sedTIMESTEP/%f/g\'  \
        -e \'s/sedSTEPS/%d/g\'  \
        -e \'s/sedPRINT/%d/g\'  \
        -e \'s/sedNATOMS/%d/g\'  \
        -e \'s/sedLBOX/%f %f %f/g\'  \
        -e \'s/sedRCUT/%f/g\'  \
        -e \'s/sedFORCEFIELD/%s_FF.inc/g\'  \
        -e \'s/sedPSF/%s_PSF.inc/g\'  \
        -e \'s/sedCONSTRAINT/%s_CONSTRAINT.inc/g\'  \
        FIST_TEMPLATE > run.inp' \
        % ( timestep, nsteps, print_frq, natoms, sizebox[0], sizebox[1], sizebox[2], \
        rcut, molecule, molecule, molecule) 
    os.system(cmd)


def psf( size3D, initcharge, mol_name):
    fileout = open(mol_name + '_PSF.inc', 'w')
    number_mol = float(size3D[0])*float(size3D[1])*float(size3D[2])
    pos_mol    = (float(initcharge[0]) - 1 )*float(size3D[1])*float(size3D[2]) + \
                 (float(initcharge[1]) - 1 )*float(size3D[2]) + \
                 float(initcharge[2]) 
    for mol in range(int(number_mol)):
        if ( (mol+1) == pos_mol):
           name = mol_name + '_CHARGE.psf'
        else:
           name = mol_name + '_NEUTRE.psf'
        result="""                               
            &MOLECULE
               NMOL              1
               CONN_FILE_NAME    ./%s
               CONN_FILE_FORMAT  PSF
            &END MOLECULE""" % name
        fileout.write(result)

def constraint(nmol, natom_mol, mol_name):
    fileout = open(mol_name + '_CONSTRAINT.inc', 'w')
    list_fix = ' '
    for imol in range(nmol):
        index = 1 + imol*natom_mol
        list_fix = list_fix + '   ' + str(index)
    result="""
            &FIXED_ATOMS
               COMPONENTS_TO_FIX XYZ
               LIST    %s
             &END FIXED_ATOMS    
    """ % list_fix
    fileout.write(result)

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



def name_bucket(BUCKET, is_a_test):
    if not is_a_test:
       complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
       md5 = hashlib.md5()
       md5.update(complete_time)
       md5name = md5.hexdigest()
    else:
       md5name = 'TEST'

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
    return bucket_name

def mkdir_bucket(bucket_name):
    test = os.path.exists(bucket_name)
    if (test):
       print "THE DIRECTORY EXISTS ALREADY!"
       sys.exit()
    else:
       os.mkdir(bucket_name)

def mkdir(file_name):
    test = os.path.exists(file_name)
    if (not test):
       os.mkdir(file_name)
    
def checkdir(file_name):
    test = os.path.exists(file_name)
    if (not test):
       print "THE DIRECTORY %s DOES NOT EXIST!" % file_name
       sys.exit()


MODULEPATH = '/scratch/grudorff/modulefiles module load mpich'
PWD = os.getcwd()
PATH_EXE = '/scratch/grudorff/antoine/bin'
PATH_TEMPLATE = PWD  + '/TEMPLATE/'
PATH_CONFIG = PWD + '/CONFIG_INITIAL'
exe_cp2k = PATH_EXE + '/cp2k.sopt'











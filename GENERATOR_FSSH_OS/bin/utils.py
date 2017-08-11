import string, re, struct, sys, math, os, time
import hashlib
import subprocess
from numpy import dot, array, prod, power
from numpy.linalg import norm
import numpy as np
from shutil import copyfile

sed_dict = {'ENSEMBLE': 'NVE',
            'TIMESTEP': 0.5,
            'TEMPERATURE': 298,
            'PRINT': 5,
            'ALPHA': 17,
            'OSPLINE': 6,
            'GMAX': 17,
            'LBOX': 100,
            'RESTRAINT': 0.005,
            'NORBITALS': 1,
            'CUTOFF_SITES': 12.0,
            'CUTOFF_CONN': 3.50,
            'SCALING': 0.0195,
            'DECO_CRIT': 1E-06,
            'CBAR': 0.50820,
            'CUTOFF_OVERLAP': 1.0E-17,
            'ELECTRONIC_STEPS': 5,
            'ANALYTICS': 'F',
            'METHOD_RESCALING': 'NACV',
            'METHOD_ADIAB_NACV': 'FAST',
            'METHOD_REVERSAL': 'NEVER',
            'NACV_INCREMENT': 1.8872589E-3,
            'PROPAGATION': 'FSSH',
            'PERIODIC': 'NONE',
            'CENTER_OF_MASS': 'T',
            'SELECT_FIRST_ADIABAT' : 'F',
            'SURF_HOP_CHOICE' : 'TRIVIAL_HOP_CORRECT',
            'T_THRESHOLD' : 0.003,
            'DECOHERENCE_CORRECTIONS': 'DAMPING',
            'EDC_C': 1.0,
            'EDC_E0': 0.1,
            'FIRST_ADIABAT': 1,
            'FIRST_DIABAT': 1,
            'INITIALIZATION': 'DIABATIC',
            'CC_CHARGED' :   1.4008,
            'CONSTRAINT_LENGTH' : 6,
            'PRINT_FSSH'        : 1,
            'K_CC_CHARGED'      : 0.263099,
	        'SEED'              : 2000,
            'REPRESENTATION'    : 'DIABATIC_BASIS',
            'RK_PROPAGATION'    : 'DIABATIC_RK'
            }



def abc_to_hmatrix(a, b, c, alpha, beta, gamma):
    """ Box vectors from box vector lengths and box vector angles (in degrees)."""
    alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
    result = np.zeros((3, 3))

    a = np.array((a, 0, 0))
    print a
    print alpha, beta, gamma
    b = b * np.array((math.cos(gamma), math.sin(gamma), 0))
    bracket = (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    c = c * np.array((math.cos(beta), bracket, math.sqrt(math.sin(beta) ** 2 - bracket ** 2)))

    result[:, 0] = a
    result[:, 1] = b
    result[:, 2] = c
    print result
    return a, b, c




class Dir(object):
    """
    """

    def __init__(self, name, paths={}, target = None):
        self._name = name
        self.path = os.getcwd() + '/' + self._name + '/'
        if target is None:
            paths.update({self._name: self.path})
        else:
            paths.update({ target : self.path})

    def clean(self):
        test = os.path.exists(self.path)
        if (test):
            os.system('rm %s/*.tmp ' % self.path)

    def checkdir(self):
        test = os.path.exists(self.path)
        if (not test):
            print "THE DIRECTORY %s DOES NOT EXIST: %s!" % (self._name, self.path)
            sys.exit()

    def mkdir(self):
        test = os.path.exists(self.path)
        if (not test):
            os.mkdir(self.path)

    def rm_mkdir(self):
        test = os.path.exists(self.path)
        if (not test):
            os.mkdir(self.path)
        else:
            os.system('rm -r %s' % self.path)
            os.mkdir(self.path)
            print " %s EXISTED, WAS REMOVED AND CREATED AGAIN" % self._name

    def chdir(self):
        self.checkdir()
        os.chdir(self.path)


class Bucket(object):
    """
        """

    def __init__(self, bucket_dict):
        self._method = bucket_dict.get('KIND_RUN')
        if bucket_dict.get('MOL_NAME') is not None:
            self._name = self._method + '-' + bucket_dict.get('MOL_NAME')
        elif bucket_dict.get('TITLE') is not None:
            self._name = self._method + '-' + bucket_dict.get('TITLE')
        else:
            self._name = self._method
        is_test = bucket_dict.get('TEST')
        if is_test == 'YES':
            self.is_test = True
        elif is_test == 'NO':
            self.is_test = False
        else:
            print "TEST must be YES or NO."
            sys.exit()

    def _get_md5name(self):
        complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
        md5 = hashlib.md5()
        md5.update(complete_time)
        md5name = md5.hexdigest()
        return md5name

    def name(self):
        if not self.is_test:
            md5name = self._get_md5name()
            # else:
            #   md5name = 'TEST'
            short_time = time.strftime("%y%m%d", time.localtime())
            bucket_name = self._name + '-' + short_time + '-' + md5name
            old_name = os.path.basename(os.getcwd())
            if (bucket_name != old_name):
                os.chdir('..')
                os.system('mv %s %s' % (old_name, bucket_name))
                os.chdir(bucket_name)
        self.path = os.getcwd()

    def checkdir(self, adir):
        test = os.path.exists(adir)
        if (not test):
            print "THE DIRECTORY %s DOES NOT EXIST!" % adir
            sys.exit()

    def mkdir(self, adir):
        test = os.path.exists(adir)
        if (not test):
            os.mkdir(adir)


class InputFile(object):
    """
    """

    def __init__(self, name):
        self._name = name
        if self._name == 'NONE':
            self.dict = {}
        else:
            self._input = open(name, 'r')
            self._lines = self._input.readlines()
            self._input.close()
            self._read()

    def _read(self):
        self.dict = {}
        for line in self._lines:
            l = string.strip(line.rstrip('\r\n'))
            info = re.split('\s+', l)
            key = info.pop(0)
            value = []
            for val in info:
                if self._convert_float(val) is None:
                    value.append(val)
                else:
                    if self._convert_int(val) is None:
                        value.append(self._convert_float(val))
                    else:
                        value.append(self._convert_int(val))
            if len(value) == 1:
                self.dict.update({key: value.pop(0)})
            else:
                self.dict.update({key: value})

    def _convert_float(self, s):
        try:
            return float(s)
        except ValueError:
            return None

    def _convert_int(self, s):
        try:
            return int(s)
        except ValueError:
            return None


class OSCrystal(object):
    """
        """

    def __init__(self, structure_dict, paths):
        self.paths = paths
        self._filemol = structure_dict.get('FILE_UNIT')
        self._sizecrystal = structure_dict.get('SIZE_CRYSTAL')
        #self._coordcharge = structure_dict.get('COORD_CHARGE')
        self._filecrystal = structure_dict.get('FILE_CRYSTAL')
        self._vecta, self._vectb, self._vectc = abc_to_hmatrix( structure_dict['ABC'][0], structure_dict['ABC'][1], structure_dict['ABC'][2],
                        structure_dict['ALPHA_BETA_GAMMA'][0], structure_dict['ALPHA_BETA_GAMMA'][1],
                        structure_dict['ALPHA_BETA_GAMMA'][2]
                        )
        print structure_dict['ABC'][0], structure_dict['ABC'][1], structure_dict['ABC'][2],\
                        structure_dict['ALPHA_BETA_GAMMA'][0], structure_dict['ALPHA_BETA_GAMMA'][1],\
                        structure_dict['ALPHA_BETA_GAMMA'][2]
        print self._vecta, self._vectb, self._vectc



    def construct_organic_crystal(self):
        with open('%s/%s' % (self.paths['structures'], self._filemol), 'r') as molfile:
            if '.xyz' in self._filemol:
                lines = molfile.readlines()[2:]
            else:
                lines = molfile.readlines()
        self._number_atoms_unit = len(lines)
        with open('%s/%s' % (self.paths['crystal'], self._filecrystal), 'w') as crystalfile:
            if '.xyz' in self._filecrystal:
                header = """%s
blabla
""" % self._calculate_number_atoms()
                crystalfile.write(header)
            for mol0_index in range(self._sizecrystal[0]):
                for mol1_index in range(self._sizecrystal[1]):
                    for mol2_index in range(self._sizecrystal[2]):
                        index3d = [mol0_index + 1, mol1_index + 1, mol2_index + 1]
                        for line_index in range(len(lines)):
                            l = string.strip(lines[line_index])
                            info = re.split('\s+', l)
                            atom_label = info[0]
                            atom_coord = [float(info[1]), float(info[2]), float(info[3])]
                            vec_shift = mol0_index * array(self._vecta) + \
                                        mol1_index * array(self._vectb) + \
                                        mol2_index * array(self._vectc)
                            result = '%s  %s\n' \
                                     % (atom_label, str(atom_coord + vec_shift).strip('[]'))
                            crystalfile.write(result)



    def _calculate_number_atoms(self):
        return self._number_atoms_unit * self._sizecrystal[0] * self._sizecrystal[1] * self._sizecrystal[2]


class CP2KRun(object):
    """
    """

    def __init__(self, dict, paths, **kwargs):
        self.paths = paths
        self._my_sed_dict = sed_dict
        self._my_sed_dict.update(dict)
        self._my_sed_dict.update(**kwargs)
        self._template_file = self._my_sed_dict.get('TEMPLATE_FILE')
        self._timestep = self._my_sed_dict.get('TIMESTEP')
	self._parallel = self._my_sed_dict.get('PARALLEL', False)

    def print_info(self):
        pass

    def _get_templates(self):
        pass

    def _get_coord(self):
        if self._my_sed_dict.get('RESTART') is not None:
            self._use_restart(self._my_sed_dict.get('RESTART'))
        else:
            self._get_new_coord()

    def _get_new_coord(self):
        pass

    def _use_restart(self, ndir):
        pass

    def _clean_velocities(self, fileinname, fileoutname):
        filein = open(fileinname)
        fileout = open(fileoutname, 'w')
        for line in filein.readlines():
            fileout.write( '\t'.join(map(str, line.split()[1:])) + '\n')
        filein.close()
        fileout.close()

    def _create_velocities(self, filevelname, filecoordname):
        filevel = open(filevelname, 'w')
        filecoord = open(filecoordname)
        for line in filecoord.readlines():
            filevel.write(' 0.000   0.000  0.000\n')
        filecoord.close()
        filevel.close()

    def run(self, ndir):
        self.ndir = ndir
        dir = Dir('run-%d' % ndir)
        self._dir = dir
        dir.rm_mkdir()
        self._complete_dict()
        self._write_input()
        self.print_info()
        self._get_input(dir)
        complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
        print "CP2K STARTS AT: " + complete_time
        dir.chdir()
	if self._parallel:
        	val = os.system(self.paths.get('cp2k') + '  run.inp > run.log')
	else:
        	#inputfile = open('run.inp')
        	logfile = open('run.log', 'w')
        	val = subprocess.call([self.paths.get('cp2k'), '-i', 'run.inp'], stdout = logfile)
        	#val = subprocess.Popen([ self.paths.get('cp2k'), '-i', 'run.inp', '-o', 'run.log'  ])

        	#val = subprocess.call([self.paths.get('cp2k')], stdin=input, stdout=logfile)
        os.chdir(self.paths.get('bucket'))
        complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
        print "CP2K FINISHES AT: " + complete_time
        if val != 0:
            fail = Dir('fail-%d' % ndir)
            fail.rm_mkdir()
            os.system(' mv %s/* %s ' % (dir.path, fail.path))
            os.rmdir(dir.path)
        return ndir + 1

    def _get_input(self, dir):
        os.system('mv %srun.inp %s' % (self.tmp.path, dir.path))
        os.system('cp %s/*psf %s' % (self.tmp.path, dir.path))

    def _write_input(self):
        self.tmp = Dir('tmp-%d' % self.ndir)
        self.tmp.rm_mkdir()
        self._get_templates()
        self._get_coord()
        os.chdir(self.tmp.path)
        self._write_topo()
        os.chdir(self.paths.get('bucket'))

    def _write_file(self, namein, nameout, number=1):
        filein = open(namein)
        fileout = open(nameout, 'w')
        result = self._amend(filein, number)
        fileout.write(result)
        filein.close()
        fileout.close()

    def _write_topo(self):
        self._write_file(self._template_file, 'run.inp')

    def _complete_dict(self):
        pass

    def _amend(self, file, number=1):
        result = "\n"
        for mol in range(1, number + 1):
            file.seek(0)
            for line in file.readlines():
                if ('INCLUDE' in line) and not ('@' in line):
                    result = result + self._include(line, mol)
                elif 'sed' in line:
                    result = result + self._sed(line)
                else:
                    result = result + line
        return result

    def _include(self, line, mol=0):
        word_list = line.split()
        if word_list[0] != 'INCLUDE':
            print 'PROBLEM IN INCLUDE FUNCTION'
            sys.exit()
        if word_list[1] == 'SPECIAL':
            file_to_include = open(word_list[2] + '-' + str(mol) + '.tmp')
        else:
            file_to_include = open(word_list[1])
        return self._amend(file_to_include)
        #result = ''
        #for my_line in file_to_include.readlines():
        #    result = result + my_line
        #return result

    def _sed(self, line):
        word_list = line.split()
        result = line
        for word in word_list:
            line = result
            if word[0:3] == 'sed':
                key = word[3:len(word)]
                val = self._my_sed_dict.get(key, '!!!! NO DEFAULT VALUE !!!')
                result = line.replace(word, str(val), 1)
        return result


class FSSHOSCrystal(CP2KRun):
    """
    """

    def __init__(self, dict, paths, **kwargs):
        super(FSSHOSCrystal, self).__init__(dict, paths, **kwargs)
        self._init = self._my_sed_dict.get('INIT')
        self._printfrq = self._my_sed_dict.get('PRINTFRQ')
        self._sizecrystal = self._my_sed_dict.get('SIZE_CRYSTAL')
        self._coordcharge = self._my_sed_dict.get('COORD_CHARGE')
        self._built_list_activated()
        self._mol_name = self._my_sed_dict.get('MOL_NAME')
        self._template_file = self._my_sed_dict.get('TEMPLATE_FILE')
        self._forcefield_file = '%s_FF.inc'  % (self._mol_name)
        self._forcefield_format = '.inc'
        if self._my_sed_dict.get('FORCEFIELD_FILE'):
            self._forcefield_file = self._my_sed_dict.get('FORCEFIELD_FILE')
            if 'prm' in self._forcefield_file:
                self._forcefield_format = '.prm'
        print self._forcefield_file
        print self._my_sed_dict.get('FORCEFIELD_FILE')
        self._filemol = 'COORD.tmp'
        self._restraint = self._my_sed_dict.get('RESTRAINT')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        self._initial_path = paths.get('initial')
        self._natom_mol = self._my_sed_dict.get('NATOM_MOL')
        print self._my_sed_dict.get('VELOCITIES')
        if self._my_sed_dict.get('VELOCITIES') is not None:
            self._do_velocities = self._my_sed_dict.get('VELOCITIES')
        else:
            self._do_velocities = True


    def _use_restart(self, restart_dict):
        self._gather_vel_coord(restart_dict, self._dir.path)
        self._clean_velocities('%s/vel-init.xyz' % self._dir.path, '%s/VELOC.init' % self._dir.path)
        #os.system('tail -%d run-%d/run-pos-1.xyz > %s/pos-init.xyz' % (self._my_sed_dict.get('NATOMS') + 2, ndir, self._dir.path) )
        #os.system('tail -%d run-%d/run-vel-1.xyz > %s/VELOC.tmp' % (self._my_sed_dict.get('NATOMS'), ndir, self._dir.path) )
        #self._clean_velocities('%s/VELOC.tmp' % self._dir.path,'%s/VELOC.init' % self._dir.path)

    def _clean_velocities(self, fileinname, fileoutname):
        filein = open(fileinname)
        fileout = open(fileoutname, 'w')
        for line in filein.readlines()[2:]:
            fileout.write( '\t'.join(map(str, line.split()[1:])) + '\n')
        filein.close()
        fileout.close()


    def _check_line(self, string, line):
        #print string
        if line == -1:
            return True
        else:
            if (re.findall(r"i = *([0-9]*)", string)):
                if  int(re.findall(r"i = *([0-9]*)", string)[0]) == line:
                    print 'Use config: %s' % line
                    return True
                else:
                    return False

    def _check_activate_mol(self, atom):
        molecule = int( atom / self._my_sed_dict['NATOM_MOL']) + 1 # WARNING: _list_activated in 1-basis
        if molecule in self._list_activated:
            return True
        else:
            return False

    def _modified(self, line):
        atom, x, y, z = line.split()
        if atom == 'C':
            label = 'CP'
        else:
            label = atom
        return '%s   %s  %s  %s\n'  % (label, x, y, z)


    def _gather_vel_coord(self, restart_dict, path):
        for iprop in ['vel', 'pos']:
            with open('%s/run-%s-1.xyz'  % (restart_dict['RESTART_DIR'], iprop), 'r') as file_:
                activate = False
                atom = 0
                with open('%s/%s-init.xyz' % (path, iprop), 'w') as fileout:
                    fileout.write('%s\n' % self._my_sed_dict['NATOMS'])
                    for ind, line in enumerate(file_.readlines()):
                        if self._check_line(line, restart_dict['CONFIG']):
                            fileout.write(line)
                            activate = True
                        elif activate:
                            if self._check_activate_mol(atom):
                                fileout.write(self._modified(line))
                            else:
                                fileout.write(line)
                            atom = atom + 1
                            if atom > (self._my_sed_dict['NATOMS'] - 1):
                                break

    def _built_list_activated(self):
        self._list_activated = self._my_sed_dict['LIST_ACTIVATED']

    def _get_input(self, dir):
        os.system('cp %s/*.txt %s' % (self.paths.get('topologies'), dir.path))

    def _write_input(self):
        self._get_coord()
        self._write_topo()

    def _write_topo(self):
        #self._write_file('%s/%s' % (self.paths['topologies'], self._forcefield_file), '%s/FORCEFIELD.include' % self._dir.path) # CREATE FORCEFIELD.include
        self._force_eval() # CREATE FORCE_EVAL.include
        self._topology() # CREATE TOPOLOGY,include
        self._aom()
        self._complete_main_input()

    def _complete_main_input(self):
        self._write_file('%s/%s' % (self.paths['templates'], self._template_file), '%s/run.inp' % self._dir.path)
        with open('%s/run.inp' % self._dir.path, 'ab+') as file_:
            for molecule in self._list_activated:
                result = "@SET  ACTIVE_MOL %s\n" % molecule
                result += "@INCLUDE FORCE_EVAL.include\n"
                file_.write(result)
        #self._write_file(self._tem, 'FORCEEVAL.tmp', number=len(self._list_activated))

    def _kind(self):
        result = """
                &KIND CP
                        ELEMENT C
                &END KIND
                &KIND H
                        ELEMENT H
                &END KIND
                &KIND CN
                        ELEMENT C
                &END KIND
        """
        return result

    def _new_psf(self):
        result = ""
        index = 0
        for molecule in self._list_activated:
            if (molecule - index - 1) != 0:
                result += """
                            &MOLECULE
                                NMOL              %s
                                CONN_FILE_NAME    ../topologies/%s
                                CONN_FILE_FORMAT  UPSF
                           &END MOLECULE
                   """ % \
                      ( molecule - index - 1, self._mol_name + "_NEUTRE.psf")
            result += """
                           @IF ${ACTIVE_MOL} == %s
                                &MOLECULE
                                    NMOL              1
                                    CONN_FILE_NAME    ../topologies/%s
                                    CONN_FILE_FORMAT  UPSF
                               &END MOLECULE
                            @ENDIF
                            @IF ${ACTIVE_MOL} /= %s
                                &MOLECULE
                                    NMOL              1
                                    CONN_FILE_NAME    ../topologies/%s
                                    CONN_FILE_FORMAT  UPSF
                               &END MOLECULE
                            @ENDIF
                   """ % \
                      ( molecule, self._mol_name + "_CHARGE.psf",
                        molecule, self._mol_name + "_NEUTRE.psf")
            index = molecule
        if (index != self._nmol):
            result += """
                            &MOLECULE
                                NMOL              %s
                                CONN_FILE_NAME    ../topologies/%s
                                CONN_FILE_FORMAT  UPSF
                            &END MOLECULE
                    """ % \
                      ( self._nmol - index, self._mol_name + "_NEUTRE.psf")
        return result

    def _amend_text(self, text, number=1):
        result = ""
        for mol in range(1, number + 1):
            for line in text.split('\n'):
               # if ('INCLUDE' in line) and not ('@' in line):
               #     result = result + self._include(line, mol)
                if 'sed' in line:
                    result = result + self._sed(line)
                else:
                    result = result + line
                result += "\n"
        return result

    def _topology(self):
        self._kind()
        with open('%s/TOPOLOGY.include' % self._dir.path, 'w') as file_:
            result = """
                &CELL
                        ABC                    %s
                        ALPHA_BETA_GAMMA       %s
                        MULTIPLE_UNIT_CELL     %s
                        PERIODIC               XYZ
                &END CELL
                %s
                &TOPOLOGY
                        COORD_FILE_FORMAT XYZ
                        COORD_FILE_NAME   pos-init.xyz
                        NUMBER_OF_ATOMS                   sedNATOMS
                        CONN_FILE_FORMAT  MOL_SET
                        &MOL_SET
                %s
                        &END MOL_SET
                        &DUMP_PSF
                                &EACH
                                   MD          1
                                &END EACH
                                FILENAME ./input
                        &END DUMP_PSF
                &END TOPOLOGY
            """ %\
            (
             '    '.join(map(str, self._my_sed_dict['ABC'])),
             '    '.join(map(str, self._my_sed_dict['ALPHA_BETA_GAMMA'])),
             '    '.join(map(str, self._my_sed_dict['SIZE_CRYSTAL'])),
             self._kind(),
             self._new_psf() )
            file_.write( self._amend_text(result))

    def _get_new_coord(self):
        os.system('cp %s/pos-%d.init %s/COORD.init' % (self._initial_path, self._init, self._dir.path))
        self._clean_velocities('%s/vel-%d.init' % (self._initial_path, self._init), '%s/VELOC.init' % self._dir.path)

    def _complete_dict(self):
        self._my_sed_dict.update({
            'NMOL': prod(self._my_sed_dict.get('SIZE_CRYSTAL'))*self._my_sed_dict['NMOL_UNIT'],
            'NDIABAT': (len(self._list_activated) ) * self._my_sed_dict.get('NORBITALS')
        })
        if self._my_sed_dict.get('RCUT') is None:
            self._my_sed_dict.update({'RCUT': 12})
        self._my_sed_dict.update({
            'FORCE_EVAL_ORDER' : '1..%d' % (len(self._list_activated) + 1),
            'NATOMS'           : self._my_sed_dict['NMOL']*self._my_sed_dict['NATOM_MOL']
        })
        self._nmol = self._my_sed_dict.get('NMOL')


    def _aom(self):
        for mol in range(1, 1 +self._nmol):
            if mol in self._list_activated:
                os.system('cat %s/%s_AOM.inc >> %s/AOM_COEFF.include' % (self.paths.get('topologies'), self._mol_name, self._dir.path))
            else:
                with open('%s/AOM_COEFF.include' % (self._dir.path), 'ab+') as file_:
                    for atom in range(self._natom_mol):
                        file_.write('XX   1    0   0.0   0.0\n')

    def _force_eval(self):
        if self._forcefield_format == '.prm':
            at_include = """
            PARMTYPE CHM
			PARM_FILE_NAME ../topologies/%s
            """ % self._forcefield_file
        else:
            at_include = """
            @INCLUDE FORCEFIELD.include
            """
        if self._do_velocities:
            velocity = """
        &VELOCITY
            @INCLUDE VELOC.init
        &END VELOCITY
"""
        else:
            velocity = "# No initial velocities"

        result = """
&FORCE_EVAL
    METHOD  FIST
    &MM
        &FORCEFIELD
            &SPLINE
                RCUT_NB         sedRCUT
            &END SPLINE
%s
        &END FORCEFIELD
        &POISSON
            &EWALD
                #EWALD_TYPE  NONE
                EWALD_TYPE  EWALD
                ALPHA        sedALPHA
                GMAX         sedGMAX
                O_SPLINE     sedOSPLINE
            &END EWALD
        &END POISSON
        &PRINT
            &ITER_INFO  SILENT
            &END ITER_INFO
            &PROGRAM_RUN_INFO LOW
                &EACH
                     MD 1
                &END EACH
            &END PROGRAM_RUN_INFO
        &END PRINT
    &END MM
    &SUBSYS
        #&COORD
        #    @INCLUDE COORD.init
        #&END COORD
        %s
        @INCLUDE TOPOLOGY.include
    &END SUBSYS
&END FORCE_EVAL
            """ % (at_include, velocity)
        with open('%s/FORCE_EVAL.include' % self._dir.path, 'w') as file_:
                file_.write( self._amend_text(result))



class FISTOSCrystal(FSSHOSCrystal):
    """
    """

    def __init__(self, dict, paths, **kwargs):
        super(FISTOSCrystal, self).__init__(dict, paths, **kwargs)


    def _get_new_coord(self):
        #os.system('tail -%s %s/crystal.xyz > %s/COORD.init' % (self._my_sed_dict['NATOMS'], self.paths['crystal'], self._dir.path))
        os.system('cp %s/crystal.xyz  %s/pos-init.xyz' % (self.paths['crystal'], self._dir.path))


    def _built_list_activated(self):
        self._list_activated = [ (int((float(self._coordcharge[0]) - 1) * float(self._sizecrystal[1] * self._sizecrystal[2]) + \
               (float(self._coordcharge[1]) - 1) * float(self._sizecrystal[2]) + \
               float(self._coordcharge[2])) - 1)*self._my_sed_dict['NMOL_UNIT'] + 1 + 1]
        print "Number charged molecule: ", self._list_activated

    def _aom(self):
        pass


class FISTOSNeutralCrystal(FISTOSCrystal):
    """
    """

    def __init__(self, dict, paths, **kwargs):
        super(FISTOSNeutralCrystal, self).__init__(dict, paths, **kwargs)


    def _built_list_activated(self):
        self._list_activated = []


    def _complete_main_input(self):
        self._write_file('%s/%s' % (self.paths['templates'], self._template_file), '%s/run.inp' % self._dir.path)
        with open('%s/run.inp' % self._dir.path, 'ab+') as file_:
            result = "@INCLUDE FORCE_EVAL.include\n"
            file_.write(result)

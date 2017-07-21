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
            'TRIVIAL' : 'TRIVIAL_HOP_CORRECT',
            'T_THRESHOLD' : 0.003,
            'DECO': 'NO_DECO_CORR',
            'EDC_C': 1.0,
            'EDC_E0': 0.1,
            'FIRST_ADIABAT': 1,
            'FIRST_DIABAT': 1,
            'INITIALIZATION': 'DIABATIC',
            'CC_CHARGED' :   1.4008,
            'CONSTRAINT_LENGTH' : 6,
            'PRINT_FSSH'        : 1,
            'K_CC_CHARGED'      : 0.263099,
	        'SEED'              : 2000
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
        self._coordcharge = structure_dict.get('COORD_CHARGE')
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
        #val = os.system(self.paths.get('cp2k') + '  run.inp > run.log')
        inputfile = open('run.inp')
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


    def _built_list_activated(self):
        _vector_diff = np.array(self._my_sed_dict.get('FIRST_MOL_CHAIN')) \
                       - np.array(self._my_sed_dict.get('LAST_MOL_CHAIN'))
        n = 0
        self._list_activated = []
        for mol0_index in range(self._sizecrystal[0]):
            for mol1_index in range(self._sizecrystal[1]):
                for mol2_index in range(self._sizecrystal[2]):
                    n += 1
                    index3d = [mol0_index + 1, mol1_index + 1, mol2_index + 1]
                    index3diff = np.array(index3d) \
                                 - np.array(self._my_sed_dict.get('FIRST_MOL_CHAIN'))
                    index3difflast = np.array(index3d) \
                                 - np.array(self._my_sed_dict.get('LAST_MOL_CHAIN'))
                    if ((np.cross(index3diff, _vector_diff)) == np.array([0, 0, 0]) ).all():
                        self._list_activated.append(n)
        print self._list_activated

    def _get_input(self, dir):
        os.system('cp %s/*.psf %s' % (self.paths.get('topologies'), dir.path))

    def _write_input(self):
        self._get_coord()
        self._write_topo()

    def _write_topo(self):
        self._write_file('%s/%s' % (self.paths['topologies'], self._forcefield_file), '%s/FORCEFIELD.include' % self._dir.path) # CREATE FORCEFIELD.include
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
                                CONN_FILE_NAME    ./%s
                                CONN_FILE_FORMAT  PSF
                           &END MOLECULE
                   """ % \
                      ( molecule - index - 1, self._mol_name + "_NEUTRE.psf")
            result += """
                           @IF ${ACTIVE_MOL} == %s
                                &MOLECULE
                                    NMOL              1
                                    CONN_FILE_NAME    ./%s
                                    CONN_FILE_FORMAT  PSF
                               &END MOLECULE
                            @ENDIF
                            @IF ${ACTIVE_MOL} /= %s
                                &MOLECULE
                                    NMOL              1
                                    CONN_FILE_NAME    ./%s
                                    CONN_FILE_FORMAT  PSF
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
                                CONN_FILE_NAME    ./%s
                                CONN_FILE_FORMAT  PSF
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
                        PERIODIC               XYZ
                &END CELL
                %s
                &TOPOLOGY
                        COORD_FILE_FORMAT XYZ
                        COORD_FILE_NAME   coord-init.xyz
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
        os.system('cp %s/crystal.xyz  %s/coord-init.xyz' % (self.paths['crystal'], self._dir.path))


    def _built_list_activated(self):
        self._list_activated = [ int((float(self._coordcharge[0]) - 1) * float(self._sizecrystal[1] * self._sizecrystal[2]) + \
               (float(self._coordcharge[1]) - 1) * float(self._sizecrystal[2]) + \
               float(self._coordcharge[2])) ]
        print self._list_activated

    def _aom(self):
        pass






class FSSHParcel(object):
    """
    """

    def __init__(self, dict, paths):
        self._my_sed_dict = sed_dict
        self._my_sed_dict.update(dict)
        self.paths = paths
        self._timestep = self._my_sed_dict.get('TIMESTEP')
        self._lengtheq = self._my_sed_dict.get('NEQ')
        self._nprod = self._my_sed_dict.get('NPROD')
        self._sizecrystal = self._my_sed_dict.get('SIZE_CRYSTAL')
        self._coordcharge = self._my_sed_dict.get('COORD_CHARGE')
        self._mol_name = self._my_sed_dict.get('MOL_NAME')
        self._template_file = self._my_sed_dict.get('TEMPLATE_FILE')
        self._filemol = self._my_sed_dict.get('FILEMOL')
        self._system = self._my_sed_dict.get('SYSTEM')
        self._restraint = self._my_sed_dict.get('RESTRAINT')
        self._templates_path = paths.get('templates')
        self._bucket_path = paths.get('bucket')
        self._output_path = paths.get('output')
        self._bin_path = paths.get('bin')
        self._complete_dict()

    def _complete_dict(self):
        norm_a = norm(array(self._my_sed_dict.get('VECTA')))
        norm_b = norm(array(self._my_sed_dict.get('VECTB')))
        norm_c = norm(array(self._my_sed_dict.get('VECTC')))
        norm_max = max(norm_a, norm_b, norm_c)
        n_a = self._my_sed_dict.get('SIZE_CRYSTAL')[0]
        n_b = self._my_sed_dict.get('SIZE_CRYSTAL')[1]
        n_c = self._my_sed_dict.get('SIZE_CRYSTAL')[2]
        n_max = max(n_a, n_b, n_c)
        self._my_sed_dict.update({
            'NMOL': prod(self._my_sed_dict.get('SIZE_CRYSTAL')),
            'NATOM_MOL': sum(1 for line in open(self.paths.get('templates') + self._my_sed_dict.get('FILEMOL'))),
            'NORM_LATTICE': [norm_a, norm_b, norm_c],
            'FIRST_DIABAT' : int(self._get_pos_mol())
        })
        if self._my_sed_dict.get('RCUT') is None:
            self._my_sed_dict.update( { 'RCUT': 5 * norm_max * n_max })
        self._natom_mol = self._my_sed_dict.get('NATOM_MOL')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        self._printfrq = self._my_sed_dict.get('PRINT')
        self._nmol = self._my_sed_dict.get('NMOL')


    def _get_pos_mol(self):
        return (float(self._coordcharge[0]) - 1) * float(self._sizecrystal[1] * self._sizecrystal[2]) + \
               (float(self._coordcharge[1]) - 1) * float(self._sizecrystal[2]) + \
               float(self._coordcharge[2])

    def gather(self, ndir):
        self._create()
        self.gather_vel_coord(ndir)
        self.gather_templates_bin()
        if (self._system == 'SOLVENT'):
            self._prepare_input_solvent()
        elif (self._system == 'CRYSTAL'):
            self._prepare_input_crystal()
        self._prepare_task()

    def _prepare_input_crystal(self):
        self._my_sed_dict.update({
            'NATOMS': prod(self._my_sed_dict.get('SIZE_CRYSTAL')) * \
                      sum(1 for line in open(self.paths.get('templates') + self._my_sed_dict.get('FILEMOL')))
        })
        self.task = """
    task = {
        'KIND_RUN'  : 'FSSH_OS',
        'TEMPLATE_FILE' : 'FSSH_CORE.template',
        'FORCEFIELD_FILE' : 'FSSH_FF.template',
        'FILE_INIT' : 'initial',
        'NUMBER_INIT' : %d,
        'NUMBER_RANDOM' : 1,
        'SYSTEM' : 'CRYSTAL',
        'MOL_NAME' : '%s',
        'NATOMS'       : %d,
        'NATOM_MOL'    : %d,
        'VECTA'        : %s,
        'VECTB'        : %s,
        'VECTC'        : %s,
        'SIZE_CRYSTAL' : %s,
        'COORD_CHARGE' : %s,
        'STEPS'        : 1,
        'PRINTFRQ'     : 1,
        'TEST'      :   'NO',
        'FIRST_DIABAT' : %s
        }
        """ % (
            self._my_sed_dict.get('NCONFIG', '!!!'),
            self._my_sed_dict.get('MOL_NAME', '!!!'),
            self._my_sed_dict.get('NATOMS', '!!!'),
            self._my_sed_dict.get('NATOM_MOL', '!!!'),
            self._my_sed_dict.get('VECTA', '!!!'),
            self._my_sed_dict.get('VECTB', '!!!'),
            self._my_sed_dict.get('VECTC', '!!!'),
            self._my_sed_dict.get('SIZE_CRYSTAL', '!!!'),
            self._my_sed_dict.get('COORD_CHARGE', '!!!'),
            self._my_sed_dict.get('FIRST_DIABAT', '!!!')
        )

    def _prepare_input_solvent(self):
        self._filecrystal = self._my_sed_dict.get('FILECRYSTAL')
        self._structure = self._my_sed_dict.get('SYSTEM')
        coordname = self.paths.get('output') + self._structure + '/' + self._filecrystal
        self._my_sed_dict.update({
            'NATOMS': sum(1 for line in open(coordname))
        })
        self.task = """
    task = {
        'KIND_RUN'  : 'FSSH_OS',
        'TEMPLATE_FILE' : 'FSSH_CORE.template',
        'FORCEFIELD_FILE' : 'FSSH_FF.template',
        'FILE_INIT' : 'initial',
        'PERIODIC' : 'XYZ',
        'NUMBER_INIT' : %d,
        'NUMBER_RANDOM' : 1,
        'SYSTEM' : 'SOLVENT',
        'MOL_NAME' :     '%s',
        'NAME_SOLVENT' : '%s',
        'SOLVENT'      : '%s',
        'NATOMS'       : %d,
        'NATOM_MOL'    : %d,
        'SIZE_BOX'     : %s,
        'VECTA'        : %s,
        'VECTB'        : %s,
        'VECTC'        : %s,
        'SIZE_CRYSTAL' : %s,
        'COORD_CHARGE' : %s,
        'STEPS'        : 1,
        'PRINTFRQ'     : 1,
        'TEST'      :   'NO',
        'RCUT'      :    %s,
        'FIRST_DIABAT' : %s
            }
        """ % (
            self._my_sed_dict.get('NCONFIG', '!!!'),
            self._my_sed_dict.get('MOL_NAME', '!!!'),
            self._my_sed_dict.get('NAME_SOLVENT', '!!!'),
            self._my_sed_dict.get('SOLVENT'),
            self._my_sed_dict.get('NATOMS', '!!!'),
            self._my_sed_dict.get('NATOM_MOL', '!!!'),
            self._my_sed_dict.get('SIZE_BOX', '!!!'),
            self._my_sed_dict.get('VECTA', '!!!'),
            self._my_sed_dict.get('VECTB', '!!!'),
            self._my_sed_dict.get('VECTC', '!!!'),
            self._my_sed_dict.get('SIZE_CRYSTAL', '!!!'),
            self._my_sed_dict.get('COORD_CHARGE', '!!!'),
            self._my_sed_dict.get('RCUT'),
            self._my_sed_dict.get('FIRST_DIABAT', '!!!')
        )

    def _check_line(self, string, line):
        #print string
        if line == -1:
            return True
        else:
            if  int(re.findall(r"i = *([0-9]*)", string)[0]) == line:
                print 'Use config: %s' % line
                return True
            else:
                return False

    def gather_vel_coord(self, ndir, output_path = None, line=-1):
        if output_path is None:
            output_path = self.initial.path
        nsolvent = self._my_sed_dict.get('NATOMS') - (self._nmol * self._natom_mol)
        os.chdir('run-%d' % ndir)
        pos_mol = (self._coordcharge[0] - 1) * (self._sizecrystal[1] * self._sizecrystal[2]) + \
                  (self._coordcharge[1] - 1) * (self._sizecrystal[2]) + \
                  (self._coordcharge[2]) - 1
        for iprop in ['vel', 'pos']:
            filein = open('run-' + iprop + '-1.xyz', 'r')
            lines = filein.readlines()
            iconfig = 0
            index = -1
            #for istep in range(0, self._nprod + 1, self._printfrq):
            for istep in range(0, self._nprod + 1, self._printfrq):
                #print istep, self._nprod + 1, self._printfrq, range(0, self._nprod + 1, self._printfrq)
                if self._check_line(lines[index + 2], line):
                    index = index + 2
                    iconfig = iconfig + 1
                    filename = iprop + '-' + str(iconfig) + '.init'
                    fileout = open(filename, 'w')
                    for imol in range(self._nmol):
                        for iatom in range(self._natom_mol):
                            index = index + 1
                            l = string.strip(lines[index])
                            info = re.split('\s+', l)
                            atom_label = self._choose_atom_label(info[0], imol=imol, icharge=pos_mol)
                            atom_xyz = [float(info[1]), float(info[2]), float(info[3])]
                            result = '%s  %s\n' \
                                     % (atom_label, str(atom_xyz).strip('[]'))
                            fileout.write(result)
                    for atom in range(nsolvent):
                        index += 1
                        fileout.write(lines[index])
                    fileout.close()
                    os.system(' mv %s %s' % (filename, output_path))
                else:
                    index = index + 2
                    for imol in range(self._nmol):
                        for iatom in range(self._natom_mol):
                            index = index + 1
                    for atom in range(nsolvent):
                        index += 1
        os.chdir(self._bucket_path)

    def gather_templates_bin(self):
        os.system('cp %s/*  %s' % (self._templates_path, self.templates.path))
        os.system('cp %s/*  %s' % (self._bin_path, self.bin.path))

    def _create(self):
        os.chdir(self._output_path)
        self.parcel = Dir('TONAME')
        self.parcel.mkdir()
        self.parcel.chdir()
        self.bin = Dir('bin')
        self.bin.mkdir()
        self.templates = Dir('templates')
        self.templates.mkdir()
        self.taskdir = Dir('task-fssh-os')
        self.taskdir.mkdir()
        self.supinitial = Dir('initial')
        self.supinitial.mkdir()
        self.supinitial.chdir()
        self.initial = Dir(self._bucket_path.split('/')[-1])
        self.initial.mkdir()
        os.chdir(self._bucket_path)

    def _choose_atom_label(self, atom, i3d=[-1, -1, -1], i3dcharge=[-1, -1, -1], imol=-1, icharge=-1):
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

    def _prepare_task(self):
        os.chdir(self.taskdir.path)
        file = open('task.py', 'w')
        result = """
 #!/usr/bin/python

# standard modules
import string, re, struct, sys, math, os, time
import numpy

# custom modules
from utils import *



def main(inputs, paths):

    %s

    inputs.update(task)

    # SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
    bucket = Bucket(inputs)
    bucket.name()
    paths.update({'bucket': bucket.path})

    task = Dir(inputs.get('INPUT_INFO'))
    paths.update( {'task' : task.path} )

    templates = Dir('templates', paths)
    templates.checkdir()

    bin = Dir('bin', paths)
    bin.checkdir()

    initial = Dir( 'initial/' + '%s' , paths)
    initial.checkdir()
    paths.update({'initial': initial.path})


    # UPLOARD ADEQUATE MODULE
    system = inputs.get('SYSTEM')
    if system == 'CRYSTAL':
        from utils import OSCluster as Structure
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import OSwSolvent as Structure
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()


    # RUN CP2K FOR number_init * number_random FSSH RUNS
    number_init = inputs.get('NUMBER_INIT', 1)
    number_random = inputs.get('NUMBER_RANDOM', 1)
    ndir= 0

    for init in range(1, number_init + 1):
        for random in range(1, number_random + 1):
            config = Config( inputs, paths, INIT = init)
            ndir = config.run(ndir)

        """ % (self.task, self._bucket_path.split('/')[-1])
        file.write(result)
        file.close()
        os.chdir(self._bucket_path)


    def create_system_info(self, output_path = None):
        if output_path is None:
            output_path = self.initial.path
        file = open('%s/system.info' % output_path, 'w')
        result = self._get_system_info()
        file.write(result)
        file.close()

    def _get_system_info(self):
        if (self._system == 'CRYSTAL'):
            result = """
        SYSTEM        CRYSTAL
        MOL_NAME      %s
        NATOMS        %d
        NATOM_MOL     %d
        VECTA         %s
        VECTB         %s
        VECTC         %s
        SIZE_CRYSTAL  %s
        COORD_CHARGE  %s
        CC_CHARGED    %s
        RESTRAINT     %s
        RCUT          %s
        """ % (
            self._my_sed_dict.get('MOL_NAME', '!!!'),
            self._my_sed_dict.get('NATOMS', '!!!'),
            self._my_sed_dict.get('NATOM_MOL', '!!!'),
            '    '.join(map(str, self._my_sed_dict.get('VECTA', '!!!'))),
            '    '.join(map(str, self._my_sed_dict.get('VECTB', '!!!'))),
            '    '.join(map(str, self._my_sed_dict.get('VECTC', '!!!'))),
            '    '.join(map(str, self._my_sed_dict.get('SIZE_CRYSTAL', '!!!'))),
            '    '.join(map(str, self._my_sed_dict.get('COORD_CHARGE', '!!!'))),
            self._my_sed_dict['CC_CHARGED'],
            self._my_sed_dict['RESTRAINT'],
            self._my_sed_dict.get('RCUT')
            )
        elif (self._system == 'SOLVENT'):
            result = """
                PERIODIC  XYZ
                SYSTEM   SOLVENT
                MOL_NAME      %s
                NAME_SOLVENT  %s
                SOLVENT       %s
                NATOMS        %d
                NATOM_MOL     %d
                SIZE_BOX      %s
                VECTA         %s
                VECTB         %s
                VECTC         %s
                SIZE_CRYSTAL  %s
                COORD_CHARGE  %s
                RCUT          %s
                CC_CHARGED    %s
                RESTRAINT     %s
                """ % (
                self._my_sed_dict.get('MOL_NAME', '!!!'),
                self._my_sed_dict.get('NAME_SOLVENT', '!!!'),
                self._my_sed_dict.get('SOLVENT'),
                self._my_sed_dict.get('NATOMS', '!!!'),
                self._my_sed_dict.get('NATOM_MOL', '!!!'),
                '    '.join(map(str, self._my_sed_dict.get('SIZE_BOX', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('VECTA', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('VECTB', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('VECTC', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('SIZE_CRYSTAL', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('COORD_CHARGE', '!!!'))),
                self._my_sed_dict.get('RCUT'),
                self._my_sed_dict['CC_CHARGED'],
            self._my_sed_dict['RESTRAINT']
            )
        elif (self._system == 'PBC_CRYSTAL'):
            result = """
                PERIODIC  XYZ
                SYSTEM   PBC_CRYSTAL
                MOL_NAME      %s
                NATOMS        %d
                NATOM_MOL     %d
                SIZE_BOX      %s
                VECTA         %s
                VECTB         %s
                VECTC         %s
                SIZE_CRYSTAL  %s
                COORD_CHARGE  %s
                RCUT          %s
                CC_CHARGED    %s
                RESTRAINT     %s
                """ % (
                self._my_sed_dict.get('MOL_NAME', '!!!'),
                self._my_sed_dict.get('NATOMS', '!!!'),
                self._my_sed_dict.get('NATOM_MOL', '!!!'),
                '    '.join(map(str, self._my_sed_dict.get('SIZE_BOX', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('VECTA', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('VECTB', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('VECTC', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('SIZE_CRYSTAL', '!!!'))),
                '    '.join(map(str, self._my_sed_dict.get('COORD_CHARGE', '!!!'))),
                self._my_sed_dict.get('RCUT'),
                self._my_sed_dict['CC_CHARGED'],
            self._my_sed_dict['RESTRAINT']
            )
        return result



class FSSHParcelBO(FSSHParcel):
    """
    """

    def __init__(self, dict, paths, output):
        super(FSSHParcelBO, self).__init__(dict, paths)
        self._nprod = self._my_sed_dict.get('STEPS')
        self.initial = output


    def _complete_dict(self):
        norm_a = norm(array(self._my_sed_dict.get('VECTA')))
        norm_b = norm(array(self._my_sed_dict.get('VECTB')))
        norm_c = norm(array(self._my_sed_dict.get('VECTC')))
        norm_max = max(norm_a, norm_b, norm_c)
        n_a = self._my_sed_dict.get('SIZE_CRYSTAL')[0]
        n_b = self._my_sed_dict.get('SIZE_CRYSTAL')[1]
        n_c = self._my_sed_dict.get('SIZE_CRYSTAL')[2]
        n_max = max(n_a, n_b, n_c)
        self._my_sed_dict.update({
            'NMOL': prod(self._my_sed_dict.get('SIZE_CRYSTAL')),
            'NORM_LATTICE': [norm_a, norm_b, norm_c],
            'FIRST_DIABAT' : int(self._get_pos_mol())
        })
        if self._my_sed_dict.get('RCUT') is None:
            self._my_sed_dict.update( { 'RCUT': 5 * norm_max * n_max })
        self._natom_mol = self._my_sed_dict.get('NATOM_MOL')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        self._printfrq = self._my_sed_dict.get('PRINT')
        self._nmol = self._my_sed_dict.get('NMOL')





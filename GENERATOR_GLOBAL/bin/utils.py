import string, re, struct, sys, math, os, time
import hashlib
import subprocess
from numpy import dot, array, prod, power
from numpy.linalg import norm
from shutil import copyfile

sed_dict = {'ENSEMBLE': 'NVE',
            'TIMESTEP': 0.5,
            'TEMPERATURE': 298,
            'PRINT': 5,
            'ALPHA': 17,
            'OSPLINE': 6,
            'GMAX': 17,
            'LBOX': 100,
            'RESTRAINT': 1E-4,
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
            'DECO': 'NO_DECO_CORR',
            'EDC_C': 1.0,
            'EDC_E0': 0.1,
            'FIRST_ADIABAT': 1,
            'FIRST_DIABAT': 1,
            'INITIALIZATION': 'DIABATIC',
            'CC_CHARGED' :   1.4008
            }


class Dir(object):
    """
    """

    def __init__(self, name, paths={}):
        self._name = name
        self.path = os.getcwd() + '/' + self._name + '/'
        paths.update({self._name: self.path})

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


class OSCluster(object):
    """
        """

    def __init__(self, structure_dict, paths):
        self._structure = structure_dict.get('SYSTEM')
        self._filemol = structure_dict.get('FILEMOL')
        self._vecta = structure_dict.get('VECTA')
        self._vectb = structure_dict.get('VECTB')
        self._vectc = structure_dict.get('VECTC')
        self._sizecrystal = structure_dict.get('SIZE_CRYSTAL')
        self._coordcharge = structure_dict.get('COORD_CHARGE')
        self._filecrystal = structure_dict.get('FILECRYSTAL')
        self._templates_path = paths.get('templates')
        self._bucket_path = paths.get('bucket')
        self._output_path = paths.get('output')
        self._write()

    def _write(self):
        self.tmp = Dir('tmp')
        self.tmp.rm_mkdir()
        self._create_dir()
        self._print_info()
        self._get_templates()
        os.chdir(self.tmp.path)
        self._my_write()
        os.system(' mv %s %s' % (self._filecrystal, self.parcel.path))
        os.chdir(self._bucket_path)

    def _get_templates(self):
        os.system('cp %s/%s %s' % (self._templates_path, self._filemol, self.tmp.path))

    def _my_write(self):
        self._construct_organic_crystal()

    def _create_dir(self):
        os.chdir(self._output_path)
        self.parcel = Dir(self._structure)
        self.parcel.mkdir()
        os.chdir(self._bucket_path)

    def _print_info(self):
        print """\n
               MOLECULE FILE   = %s
               CRYSTAL LATTICE a = %s
               CRYSTAL LATTICE b = %s
               CRYSTAL LATTICE c = %s
               CRYSTAL SIZE    = %s
               COORD CHARGE    = %s
               CRYSTAL FILE    = %s
               \n
               """ \
              % (self._filemol,
                 str(self._vecta).strip('[]'),
                 str(self._vectb).strip('[]'),
                 str(self._vectc).strip('[]'),
                 str(self._sizecrystal).strip('[]'),
                 str(self._coordcharge).strip('[]'),
                 self._filecrystal)

    def _construct_organic_crystal(self):
        molfile = open(self._filemol, 'r')
        crystalfile = open(self._filecrystal, 'w')
        lines = molfile.readlines()
        molfile.close
        for mol0_index in range(self._sizecrystal[0]):
            for mol1_index in range(self._sizecrystal[1]):
                for mol2_index in range(self._sizecrystal[2]):
                    index3d = [mol0_index + 1, mol1_index + 1, mol2_index + 1]
                    for line_index in range(len(lines)):
                        l = string.strip(lines[line_index])
                        info = re.split('\s+', l)
                        atom_label = self._choose_atom_label(info[0], i3d=index3d, i3dcharge=self._coordcharge)
                        atom_coord = [float(info[1]), float(info[2]), float(info[3])]
                        vec_shift = mol0_index * array(self._vecta) + \
                                    mol1_index * array(self._vectb) + \
                                    mol2_index * array(self._vectc)
                        result = '%s  %s\n' \
                                 % (atom_label, str(atom_coord + vec_shift).strip('[]'))
                        crystalfile.write(result)
        crystalfile.close

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
            atom_label = None
            print "THERE IS A PROBLEM in choose_atom_label"
        return atom_label


class OSwSolvent(OSCluster):
    """
    """

    def __init__(self, inputs, paths):
        self._closest_dist = inputs.get('CLOSEST_DIST')
        self._kind_solvent = inputs.get('SOLVENT')
        self._density = inputs.get('DENSITY')
        self._sizebox = inputs.get('SIZE_BOX')
        super(OSwSolvent, self).__init__(inputs, paths)

    def _my_write(self):
        self._construct_organic_crystal()
        self._add_solvent()

    def _add_solvent(self):
        filecoord = open(self._filecrystal, 'a+')
        self._get_grid()
        self._get_carbon_pos(filecoord)
        list_of_list = [list(array(range(self._grid[i])) - int(self._grid[i] / 2)) for i in range(3)]
        grid = [[i, j, k] for i in list_of_list[0]
                for j in list_of_list[1]
                for k in list_of_list[2]]

        for pos in grid:
            realpos = array(pos) * array(self._realgrid)
            dist = self._get_os_dist(realpos)
            if dist >= self._closest_dist:
                result = "%s  %f  %f  %f \n" % (self._kind_solvent, realpos[0], realpos[1], realpos[2])
                filecoord.write(result)
        filecoord.close()

    def _get_carbon_pos(self, filecoord):
        self._carbon_pos = []
        for line in filecoord.readlines():
            if 'C' in line:
                info = re.split('\s+', line)
                self._carbon_pos.append([float(info[1]), float(info[2]), float(info[3])])

    def _get_grid(self):
        volume = prod(self._sizebox)
        pre_natoms = self._density * volume
        cubic_roots = int(power(pre_natoms, 1.0 / 3.0))
        self._grid = [cubic_roots, cubic_roots, cubic_roots]
        self._realgrid = array(self._sizebox) / array(self._grid)

    def _get_os_dist(self, pos_solvent):
        list = []
        for pos_carbon in self._carbon_pos:
            dist = norm(array(pos_carbon) - array(pos_solvent))
            list.append(dist)
        return min(list)


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
        self._complete_dict()
        self._write_input()
        self.print_info()
        dir = Dir('run-%d' % ndir)
        dir.rm_mkdir()
        os.system('mv %srun.inp %s' % (self.tmp.path, dir.path))
        os.system('cp %s/*psf %s' % (self.tmp.path, dir.path))
        complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
        print "CP2K STARTS AT: " + complete_time
        dir.chdir()
        val = os.system(self.paths.get('cp2k') + '  run.inp > run.log')
        #val = subprocess.Popen([ self.paths.get('cp2k'), '-i', 'run.inp', '-o', 'run.log'  ])
        #input = open('run.inp')
        #logfile = open('run.log', 'w')
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

    def _write_input(self):
        self.tmp = Dir('tmp-%d' % self.ndir)
        self.tmp.rm_mkdir()
        self._get_templates()
        self._get_coord()
        os.chdir(self.tmp.path)
        self._write_topo()
        self._write_file(self._template_file, 'run.inp')
        os.chdir(self.paths.get('bucket'))

    def _write_file(self, namein, nameout, number=1):
        filein = open(namein)
        fileout = open(nameout, 'w')
        result = self._amend(filein, number)
        fileout.write(result)
        filein.close()
        fileout.close()

    def _write_topo(self):
        pass

    def _complete_dict(self):
        pass

    def _amend(self, file, number=1):
        result = "\n"
        for mol in range(1, number + 1):
            file.seek(0)
            for line in file.readlines():
                if 'INCLUDE' in line:
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


class CP2KOS(CP2KRun):
    """
    """

    def __init__(self, dict, paths, **kwargs):
        super(CP2KOS, self).__init__(dict, paths, **kwargs)
        self._sizecrystal = self._my_sed_dict.get('SIZE_CRYSTAL')
        self._coordcharge = self._my_sed_dict.get('COORD_CHARGE')
        self._mol_name = self._my_sed_dict.get('MOL_NAME')
        self._name_organic = self._my_sed_dict.get('MOL_NAME')
        self._filecrystal = self._my_sed_dict.get('FILECRYSTAL')
        self._filemol = self._my_sed_dict.get('FILEMOL')
        self._structure = self._my_sed_dict.get('SYSTEM')

    def print_info(self):
        print "Hey Hey"

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
            'NATOMS': prod(self._my_sed_dict.get('SIZE_CRYSTAL')) * \
                      sum(1 for line in open(self.paths.get('templates') + self._my_sed_dict.get('FILEMOL'))),
            'NORM_LATTICE': [norm_a, norm_b, norm_c],
            'LBOXA': 10 * norm_max * n_max,
            'LBOXB': 10 * norm_max * n_max,
            'LBOXC': 10 * norm_max * n_max,
            'PRINT': int(self._my_sed_dict.get('NPROD') / self._my_sed_dict.get('NCONFIG'))
        })
        if self._my_sed_dict.get('RCUT') is None:
            self._my_sed_dict.update( { 'RCUT': 5 * norm_max * n_max })
        self._natom_mol = self._my_sed_dict.get('NATOM_MOL')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        self._restraint = self._my_sed_dict.get('RESTRAINT')

    def _use_restart(self, ndir):
        os.system('tail -%d run-%d/run-pos-1.xyz > COORD.tmp' % (self._my_sed_dict.get('NATOMS'), ndir) )
        os.system('tail -%d run-%d/run-vel-1.xyz > preVELOC.tmp' % (self._my_sed_dict.get('NATOMS'), ndir) )
        self._clean_velocities('preVELOC.tmp','VELOC.tmp')
        os.system('mv *.tmp %s' % self.tmp.path)

    def _get_new_coord(self):
        os.system('cp %s/%s/%s %s/COORD.tmp' % (self.paths.get('output'), self._structure,
                                                self._filecrystal, self.tmp.path))
        self._create_velocities('%s/VELOC.tmp'% self.tmp.path,'%s/COORD.tmp' % self.tmp.path)


    def _get_templates(self):
        os.system('cp %s/*.psf %s' % (self.paths.get('templates'), self.tmp.path))
        os.system('cp %s/*.inc %s' % (self.paths.get('templates'), self.tmp.path))
        os.system('cp %s/FIST* %s' % (self.paths.get('templates'), self.tmp.path))
        os.system('cp %s/%s %s' % (self.paths.get('templates'), self._filemol, self.tmp.path))

    def _write_topo(self):
        self._forcefield()
        self._kind()
        self._psf()
        self._colvar()
        self._constraint()

    def _forcefield(self):
        os.system('cp %s_FF.inc FORCEFIELD.tmp' % self._mol_name)

    def _kind(self):
        fileout = open('KIND.tmp', 'w')
        result = """\n
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
        fileout.write(result)
        fileout.close()

    def _constraint(self):
        fileout = open('CONSTRAINT.tmp', 'w')
        fileout.write('        &CONSTRAINT\n')
        list_mol = [[i, j, k] for i in range(1, self._sizecrystal[0] + 1)
                    for j in range(1, self._sizecrystal[1] + 1)
                    for k in range(1, self._sizecrystal[2] + 1)]
        colvar = 0
        list_mol2 = list_mol
        for mol in list_mol:
            list_mol2.remove(mol)
            if (list_mol2 is not None):
                for mol2 in list_mol2:
                    top_vect = abs(array(mol) - array(mol2))
                    top_dist = norm(top_vect)
                    if (top_dist == 1) or (top_dist == 2):
                        colvar = colvar + 1
                        result = """\n
             &COLLECTIVE
                  &RESTRAINT
                         K              %f
                  &END RESTRAINT
                  COLVAR                %d
                  INTERMOLECULAR
                  TARGET                %f
             &END COLLECTIVE
                         \n     """ \
                                 % (self._restraint,
                                    colvar,
                                    dot(array(self._norm_lattice), top_vect))
                        fileout.write(result)
        fileout.write('        &END CONSTRAINT\n')
        fileout.close()

    def _colvar(self):
        fileout = open('COLVAR.tmp', 'w')
        list_mol = [[i, j, k] for i in range(1, self._sizecrystal[0] + 1)
                    for j in range(1, self._sizecrystal[1] + 1)
                    for k in range(1, self._sizecrystal[2] + 1)]
        list_mol2 = list_mol
        for mol in list_mol:
            list_mol2.remove(mol)
            if (list_mol2 is not None):
                for mol2 in list_mol2:
                    top_vect = abs(array(mol) - array(mol2))
                    top_dist = norm(top_vect)
                    if (top_dist == 1) or (top_dist == 2):
                        result = """\n
                &COLVAR
                     &DISTANCE
                          &POINT
                               ATOMS %s
                               TYPE  GEO_CENTER
                          &END POINT
                          &POINT
                               ATOMS %s
                               TYPE  GEO_CENTER
                          &END POINT
                          POINTS 1 2
                     &END DISTANCE
                &END COLVAR
                \n
                """ % (
                            self._list_carbon(mol),
                            self._list_carbon(mol2))
                        fileout.write(result)
        fileout.close()

    def _list_carbon(self, mol):
        molfile = open(self._filemol, 'r')
        lines = molfile.readlines()
        list_carbon = []
        num_mol = (mol[0] - 1) * self._sizecrystal[1] * self._sizecrystal[2] + \
                  (mol[1] - 1) * self._sizecrystal[2] + \
                  mol[2] - 1
        for i in range(self._natom_mol):
            line = lines[i]
            if ('C' in line):
                list_carbon.append(i + 1 + (num_mol) * self._natom_mol)
        return '\t'.join(map(str, list_carbon))

    def _get_pos_mol(self):
        return (float(self._coordcharge[0]) - 1) * float(self._sizecrystal[1] * self._sizecrystal[2]) + \
               (float(self._coordcharge[1]) - 1) * float(self._sizecrystal[2]) + \
               float(self._coordcharge[2])

    def _psf(self, pos_mol=-1, count=False):
        if pos_mol == -1:
            pos_mol = self._get_pos_mol()
        if count:
            fileout = open('PSF-%d.tmp' % pos_mol, 'w')
        else:
            fileout = open('PSF.tmp', 'w')
        number_mol = prod(self._sizecrystal)
        for mol in range(int(number_mol)):
            if ((mol + 1) == pos_mol):
                name = self._name_organic + '_CHARGE.psf'
            else:
                name = self._name_organic + '_NEUTRE.psf'
            result = """\n
                    &MOLECULE
                        NMOL              1
                        CONN_FILE_NAME    ./%s
                        CONN_FILE_FORMAT  PSF
                   &END MOLECULE\n """ % name
            fileout.write(result)
        fileout.close()


class CP2KOSwSolvent(CP2KOS):
    def __init__(self, dict, paths, **kwargs):
        super(CP2KOSwSolvent, self).__init__(dict, paths, **kwargs)
        self._kind_solvent = self._my_sed_dict.get('SOLVENT')
        self._sizebox = self._my_sed_dict.get('SIZE_BOX')
        self._name_solvent = self._my_sed_dict.get('NAME_SOLVENT')

    def _complete_dict(self):
        super(CP2KOSwSolvent, self)._complete_dict()
        coordname = self.paths.get('output') + self._structure + '/' + self._filecrystal
        self._my_sed_dict.update({
            'NATOMS': sum(1 for line in open(coordname)),
            'NSOLVENT': open(coordname).read().count(self._kind_solvent),
            'LBOXA': self._sizebox[0],
            'LBOXB': self._sizebox[1],
            'LBOXC': self._sizebox[2],
            'PERIODIC': 'XYZ'
        })
        if self._my_sed_dict.get('RCUT') is None:
            self._my_sed_dict.update( { 'RCUT': min(self._sizebox) / 2 })
        self._nsolvent = self._my_sed_dict.get('NSOLVENT')

    def _kind(self):
        fileout = open('KIND.tmp', 'w')
        result = """\n
                   &KIND CP
                        ELEMENT C
                &END KIND
                &KIND H
                        ELEMENT H
                &END KIND
                &KIND CN
                        ELEMENT C
                &END KIND
                &KIND %s
                        ELEMENT %s
                &END KIND
            """ % (self._kind_solvent.title(), self._kind_solvent.title())
        fileout.write(result)
        fileout.close()

    def _psf(self, pos_mol=-1, count=False):
        if pos_mol == -1:
            pos_mol = self._get_pos_mol()
        if count:
            fileout = open('PSF-%d.tmp' % pos_mol, 'w')
        else:
            fileout = open('PSF.tmp', 'w')
        number_mol = prod(self._sizecrystal)
        for mol in range(int(number_mol)):
            if ((mol + 1) == pos_mol):
                name = self._name_organic + '_CHARGE.psf'
            else:
                name = self._name_organic + '_NEUTRE.psf'
            result = """\n
                    &MOLECULE
                        NMOL              1
                        CONN_FILE_NAME    ./%s
                        CONN_FILE_FORMAT  PSF
                   &END MOLECULE\n """ % name
            fileout.write(result)
        result = """\n
                    &MOLECULE
                        NMOL              %s
                        CONN_FILE_NAME    ./%s
                        CONN_FILE_FORMAT  PSF
                   &END MOLECULE\n """ % (self._nsolvent, self._kind_solvent.title() + '.psf')
        fileout.write(result)
        fileout.close()

    def _forcefield(self):
        print              self._mol_name + '_' + self._name_solvent
        os.system('cp %s_FF.inc FORCEFIELD.tmp' % (self._mol_name + '_' + self._name_solvent))


class CP2KOSwSolventFSSH(CP2KOSwSolvent):
    """
    """

    def __init__(self, dict, paths, **kwargs):
        super(CP2KOSwSolventFSSH, self).__init__(dict, paths, **kwargs)
        self._init = self._my_sed_dict.get('INIT')
        self._printfrq = self._my_sed_dict.get('PRINTFRQ')
        self._sizecrystal = self._my_sed_dict.get('SIZE_CRYSTAL')
        self._coordcharge = self._my_sed_dict.get('COORD_CHARGE')
        self._mol_name = self._my_sed_dict.get('MOL_NAME')
        self._template_file = self._my_sed_dict.get('TEMPLATE_FILE')
        self._forcefield_file = self._my_sed_dict.get('FORCEFIELD_FILE')
        self._filemol = 'COORD.tmp'
        self._restraint = self._my_sed_dict.get('RESTRAINT')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        self._initial_path = paths.get('initial')
        self._natom_mol = self._my_sed_dict.get('NATOM_MOL')

    def _write_topo(self):
        self._forcefield()
        self._kind()
        for mol in range(1, self._nmol + 1):
            self._psf(mol, True)
        self._colvar()
        self._constraint()
        self._aom()
        self._write_file(self._forcefield_file, 'FORCEEVAL.tmp', number=self._nmol)

    def _complete_dict(self):
        self._my_sed_dict = self._my_sed_dict
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
            'LBOXA': 10 * norm_max * n_max,
            'LBOXB': 10 * norm_max * n_max,
            'LBOXC': 10 * norm_max * n_max
        })
        self._natom_mol = self._my_sed_dict.get('NATOM_MOL')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        coordname = self.paths.get('initial') + '/pos-1.init'
        self._my_sed_dict.update({
            'NSOLVENT': open(coordname).read().count(self._kind_solvent.title()),
            'NDIABAT': prod(self._my_sed_dict.get('SIZE_CRYSTAL')) * self._my_sed_dict.get('NORBITALS'),
            'LBOXA': self._sizebox[0],
            'LBOXB': self._sizebox[1],
            'LBOXC': self._sizebox[2],
            'PERIODIC': 'XYZ'
        })
        if self._my_sed_dict.get('RCUT') is None:
            self._my_sed_dict.update( { 'RCUT': min(self._sizebox) / 2 })
        self._nsolvent = self._my_sed_dict.get('NSOLVENT')
        print "NSOLVENT", self._nsolvent
        self._nmol = self._my_sed_dict.get('NMOL')
        self._my_sed_dict.update({
            'FORCE_EVAL_ORDER': '1..%d' % (self._my_sed_dict.get('NDIABAT') + 1)
        })
        self._restraint = self._my_sed_dict.get('RESTRAINT')

    def _aom(self):
        for mol in range(self._nmol):
            os.system('cat %s_AOM.inc >> AOM_COEFF.tmp' % self._name_organic)
        file = open('AOM_COEFF.tmp', 'ab+')
        for atom in range(self._nsolvent):
            file.write('%s    1    0   0.0   0.0\n' % self._kind_solvent)
        file.close()

    def _get_new_coord(self):
        os.system('cp %s/pos-%d.init %s/COORD.tmp' % (self._initial_path, self._init, self.tmp.path))
        self._clean_velocities('%s/vel-%d.init' % (self._initial_path, self._init ),'%s/VELOC.tmp' % self.tmp.path)



    def _get_templates(self):
        os.system('cp %s/*.psf %s' % (self.paths.get('templates'), self.tmp.path))
        os.system('cp %s/*.inc %s' % (self.paths.get('templates'), self.tmp.path))
        os.system('cp %s/FSSH* %s' % (self.paths.get('templates'), self.tmp.path))


class CP2KOSFIST(CP2KOS):
    """
    """

    def __init__(self, dict, paths, **kwargs):
        super(CP2KOSFIST, self).__init__(dict, paths, **kwargs)
        self._init = self._my_sed_dict.get('INIT')
        self._printfrq = self._my_sed_dict.get('PRINTFRQ')
        self._sizecrystal = self._my_sed_dict.get('SIZE_CRYSTAL')
        self._coordcharge = self._my_sed_dict.get('COORD_CHARGE')
        self._mol_name = self._my_sed_dict.get('MOL_NAME')
        self._template_file = self._my_sed_dict.get('TEMPLATE_FILE')
        self._forcefield_file = self._my_sed_dict.get('FORCEFIELD_FILE')
        self._filemol = 'COORD.tmp'
        self._restraint = self._my_sed_dict.get('RESTRAINT')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        self._initial_path = paths.get('initial')
        self._natom_mol = self._my_sed_dict.get('NATOM_MOL')

    def print_info(self):
        print "Hey"

    def _write_topo(self):
        self._forcefield()
        self._kind()
        self._psf()
        self._colvar()
        self._constraint()


    def _get_new_coord(self):
        os.system('cp %s/pos-%d.init %s/COORD.tmp' % (self._initial_path, self._init, self.tmp.path))
        self._clean_velocities('%s/vel-%d.init' % (self._initial_path, self._init ),'%s/VELOC.tmp' % self.tmp.path)

    def _get_templates(self):
        os.system('cp %s/*.psf %s' % (self.paths.get('templates'), self.tmp.path))
        os.system('cp %s/*.inc %s' % (self.paths.get('templates'), self.tmp.path))
        os.system('cp %s/FIST* %s' % (self.paths.get('templates'), self.tmp.path))


    def _complete_dict(self):
        norm_a = norm(array(self._my_sed_dict.get('VECTA')))
        norm_b = norm(array(self._my_sed_dict.get('VECTB')))
        norm_c = norm(array(self._my_sed_dict.get('VECTC')))
        n_a = self._my_sed_dict.get('SIZE_CRYSTAL')[0]
        n_b = self._my_sed_dict.get('SIZE_CRYSTAL')[1]
        n_c = self._my_sed_dict.get('SIZE_CRYSTAL')[2]
        n_max = max(n_a, n_b, n_c)
        norm_max = max(norm_a, norm_b, norm_c)
        self._my_sed_dict.update({
            'NMOL': prod(self._my_sed_dict.get('SIZE_CRYSTAL')),
            'NDIABAT': prod(self._my_sed_dict.get('SIZE_CRYSTAL')) * self._my_sed_dict.get('NORBITALS'),
            'NORM_LATTICE': [norm_a, norm_b, norm_c],
            'LBOXA': 10 * norm_max * n_max,
            'LBOXB': 10 * norm_max * n_max,
            'LBOXC': 10 * norm_max * n_max
        })
        if self._my_sed_dict.get('RCUT') is None:
            self._my_sed_dict.update( { 'RCUT': 5 * norm_max * n_max })
        self._my_sed_dict.update({
            'FORCE_EVAL_ORDER': '1..%d' % (self._my_sed_dict.get('NDIABAT') + 1)
        })
        self._nmol = self._my_sed_dict.get('NMOL')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        self._restraint = self._my_sed_dict.get('RESTRAINT')





class CP2KOSFSSH(CP2KOS):
    """
    """

    def __init__(self, dict, paths, **kwargs):
        super(CP2KOSFSSH, self).__init__(dict, paths, **kwargs)
        self._init = self._my_sed_dict.get('INIT')
        self._printfrq = self._my_sed_dict.get('PRINTFRQ')
        self._sizecrystal = self._my_sed_dict.get('SIZE_CRYSTAL')
        self._coordcharge = self._my_sed_dict.get('COORD_CHARGE')
        self._mol_name = self._my_sed_dict.get('MOL_NAME')
        self._template_file = self._my_sed_dict.get('TEMPLATE_FILE')
        self._forcefield_file = self._my_sed_dict.get('FORCEFIELD_FILE')
        self._filemol = 'COORD.tmp'
        self._restraint = self._my_sed_dict.get('RESTRAINT')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        self._initial_path = paths.get('initial')
        self._natom_mol = self._my_sed_dict.get('NATOM_MOL')

    def print_info(self):
        print "Hey"

    def _write_topo(self):
        self._forcefield()
        self._kind()
        for mol in range(1, self._nmol + 1):
            self._psf(mol, True)
        self._colvar()
        self._constraint()
        self._aom()
        self._write_file(self._forcefield_file, 'FORCEEVAL.tmp', number=self._nmol)

    def _get_new_coord(self):
        os.system('cp %s/pos-%d.init %s/COORD.tmp' % (self._initial_path, self._init, self.tmp.path))
        self._clean_velocities('%s/vel-%d.init' % (self._initial_path, self._init ),'%s/VELOC.tmp' % self.tmp.path)

    def _get_templates(self):
        os.system('cp %s/*.psf %s' % (self.paths.get('templates'), self.tmp.path))
        os.system('cp %s/*.inc %s' % (self.paths.get('templates'), self.tmp.path))
        os.system('cp %s/FSSH* %s' % (self.paths.get('templates'), self.tmp.path))


    def _complete_dict(self):
        norm_a = norm(array(self._my_sed_dict.get('VECTA')))
        norm_b = norm(array(self._my_sed_dict.get('VECTB')))
        norm_c = norm(array(self._my_sed_dict.get('VECTC')))
        n_a = self._my_sed_dict.get('SIZE_CRYSTAL')[0]
        n_b = self._my_sed_dict.get('SIZE_CRYSTAL')[1]
        n_c = self._my_sed_dict.get('SIZE_CRYSTAL')[2]
        n_max = max(n_a, n_b, n_c)
        norm_max = max(norm_a, norm_b, norm_c)
        self._my_sed_dict.update({
            'NMOL': prod(self._my_sed_dict.get('SIZE_CRYSTAL')),
            'NDIABAT': prod(self._my_sed_dict.get('SIZE_CRYSTAL')) * self._my_sed_dict.get('NORBITALS'),
            'NORM_LATTICE': [norm_a, norm_b, norm_c],
            'LBOXA': 10 * norm_max * n_max,
            'LBOXB': 10 * norm_max * n_max,
            'LBOXC': 10 * norm_max * n_max
        })
        if self._my_sed_dict.get('RCUT') is None:
            self._my_sed_dict.update( { 'RCUT': 5 * norm_max * n_max })
        self._my_sed_dict.update({
            'FORCE_EVAL_ORDER': '1..%d' % (self._my_sed_dict.get('NDIABAT') + 1)
        })
        self._nmol = self._my_sed_dict.get('NMOL')
        self._norm_lattice = self._my_sed_dict.get('NORM_LATTICE')
        self._restraint = self._my_sed_dict.get('RESTRAINT')

    def _aom(self):
        for mol in range(self._nmol):
            os.system('cat %s_AOM.inc >> AOM_COEFF.tmp' % self._mol_name)


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
        self._local_path = paths.get('local_paths')
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

    def gather_vel_coord(self, ndir):
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
            for istep in range(0, self._nprod, self._printfrq):
                index = index + 2
                iconfig = iconfig + 1
                filename = iprop + '-' + str(iconfig) + '.init'
                fileout = open(filename, 'w')
                for imol in range(self._nmol):
                    for iatom in range(self._natom_mol):
                        index = index + 1
                        # index = (istep / self._printfrq) * (self._nmol * self._natom_mol + 2) \
                        #        + 2 \
                        #        + imol * self._natom_mol \
                        #        + iatom
                        l = string.strip(lines[index])
                        info = re.split('\s+', l)
  #                      if iprop == 'pos':
                        atom_label = self._choose_atom_label(info[0], imol=imol, icharge=pos_mol)
  #                      else:
   #                         atom_label = '   '
                        atom_xyz = [float(info[1]), float(info[2]), float(info[3])]
                        result = '%s  %s\n' \
                                 % (atom_label, str(atom_xyz).strip('[]'))
                        fileout.write(result)
                for atom in range(nsolvent):
                    # index = (self._nmol * self._natom_mol) + atom + 2  #2 from the the two initial lines at each timestep
                    index += 1
    #                if iprop == 'pos':
                    fileout.write(lines[index])
     #               else:
      #                  info = lines[index].split()
       #                 atom_xyz = [float(info[1]), float(info[2]), float(info[3])]
        #                result = ' %s\n' \
         #                        % (str(atom_xyz).strip('[]'))
          #              fileout.write(result)


                fileout.close()
                os.system(' mv %s %s' % (filename, self.initial.path))
        os.chdir(self._bucket_path)

    def gather_templates_bin(self):
        os.system('cp %s/*  %s' % (self._templates_path, self.templates.path))
        os.system('cp %s/*  %s' % (self._bin_path, self.bin.path))
        os.system('cp %s/*  %s' % (self._local_path, self.local.path))

    def _create(self):
        os.chdir(self._output_path)
        self.parcel = Dir('TONAME')
        self.parcel.mkdir()
        self.parcel.chdir()
        self.bin = Dir('bin')
        self.bin.mkdir()
        self.templates = Dir('templates')
        self.templates.mkdir()
        self.local = Dir('local_paths')
        self.local.mkdir()
        self.task = Dir('task')
        self.task.mkdir()
        self.task.chdir()
        self.subtask = Dir('fssh_os')
        self.subtask.mkdir()
        self.subtask.chdir()
        self.initial = Dir('initial')
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
        os.chdir(self.subtask.path)
        file = open('task.py', 'w')
        result = """
 #!/usr/bin/python

import string, re, struct, sys, math, os, time
import numpy

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
    templates.clean()

    bin = Dir('bin', paths)
    bin.checkdir()

    # FIND CP2K PATHS
    try:
        local_paths = Dir('local_paths', paths)
        local_paths.checkdir()
        cp2k_file = open(paths.get('local_paths') + 'cp2k.path', 'r')
        paths.update({'cp2k': cp2k_file.read().rstrip()})
        if not os.path.isfile(paths.get('cp2k')):
            raise SystemExit('WARNING: check path for CP2K executable in local_paths/cp2k.path')
    except:
        raise SystemExit("WARINING: please provide the path for CP2K executable in local_paths/cp2k.path")


    os.system(' cp -r %%s/%%s %%s' %% (paths.get('task'), inputs.get('FILE_INIT'), paths.get('bucket')))
    initial = Dir(inputs.get('FILE_INIT'), paths)
    initial.checkdir()
    paths.update({'initial': initial.path})

    system = inputs.get('SYSTEM')
    if system == 'CRYSTAL':
        from utils import OSCluster as Structure
        from utils import CP2KOSFSSH as Config
    elif system == 'SOLVENT':
        from utils import OSwSolvent as Structure
        from utils import CP2KOSwSolventFSSH as Config
    else:
        sys.exit()

    number_init = inputs.get('NUMBER_INIT', 1)
    number_random = inputs.get('NUMBER_RANDOM', 1)
    ndir= 0

    for init in range(1, number_init + 1):
        for random in range(1, number_random + 1):
            config = Config( inputs, paths, INIT = init)
            ndir = config.run(ndir)

        """ % self.task
        file.write(result)
        file.close()
        os.chdir(self._bucket_path)


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





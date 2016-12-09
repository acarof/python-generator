import string, re, struct, sys, math, os, time
import hashlib
import subprocess
from numpy import dot, array, prod
from numpy.linalg import norm
from shutil import copyfile


sed_dict = { 'ENSEMBLE'    :'NVE',
             'TIMESTEP'    : 0.5,
             'TEMPERATURE' : 298,
             'PRINT'       : 5,
             'ALPHA'       : 17,
             'OSPLINE'     : 6,
             'GMAX'        : 17,
             'RCUT'        : 50,
             'LBOX'        : 100,
             'RESTRAINT'   : 1E-4,
             'NORBITALS'   : 1,
             'CUTOFF_SITES': 12.0,
             'CUTOFF_CONN' : 3.50,
             'SCALING'     : 0.1,
             'FIRST_DIABAT': 1,
             'DECO_CRIT'   : 1E-06,
             'CBAR'        : 0.50820,
             'CUTOFF_OVERLAP' : 1.0E-17,
             'ELECTRONIC_STEPS' : 5,
             'COLLAPSE'         : 'T',
             'ANALYTICS'        : 'F',
             'METHOD_RESCALING' : 'NACV',
             'METHOD_ADIAB_NACV' : 'FAST',
             'METHOD_REVERSAL'   : 'NEVER',
             'NACV_INCREMENT'    : 1.8872589E-3,
             'PROPAGATION'       : 'FSSH'
           }

class Bucket(object):
	"""
        """
	def __init__(self, bucket_dict):
           self._name = bucket_dict.get('MOL_NAME', 'TEST')
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
            else:
               md5name = 'TEST'
            short_time = time.strftime("%y%m%d", time.localtime())
            bucket_name = self._name + '_' + short_time + '_'  '_' + md5name
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
          self.name = name
          self.input =  open( name, 'r')
          self.lines =  self.input.readlines()
          self.input.close()      
#          self.dict_list = {'BUCKET' :['NAME', 'TEST'],
#                            'CRYSTAL':['FILEMOL','VECTA','VECTB','VECTC',
#                                       'SIZE_CRYSTAL','COORD_CHARGE','FILECRYSTAL'],
#                            'SYSTEM' :['SYSTEM'],
#                            'CP2K'   :['TIMESTEP','NEQ','NPROD','PRINTFRQ','SIZE_CRYSTAL']}

      def read(self):
          self.dict = {}
          for line in self.lines:
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
              if len(value)==1: 
                 self.dict.update({key:value.pop(0)})
              else:
                 self.dict.update({key:value})

    
#      def get_info(self, info):
#          list = self.dict_list.get(info)
#          dict = { key : self._dict.get(key) for key in list}
#          return dict

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

class Crystal(object):
	"""
        """
	def __init__(self, structure_dict, paths):
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

        def print_info(self):
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
               % ( self._filemol,
               str(self._vecta).strip('[]'),
               str(self._vectb).strip('[]'),
               str(self._vectc).strip('[]'),
               str(self._sizecrystal).strip('[]'),
               str(self._coordcharge).strip('[]'),
               self._filecrystal)

        def complete_dict(self,  dict):
           if (dict) is None:
              print "Error in input.complete"
              sys.exit()
           else:
              norm_a = norm(array(dict.get('VECTA')))
              norm_b = norm(array(dict.get('VECTB')))
              norm_c = norm(array(dict.get('VECTC')))
              norm_max = max(norm_a, norm_b, norm_c)
              n_a    = dict.get('SIZE_CRYSTAL')[0]
              n_b    = dict.get('SIZE_CRYSTAL')[1]
              n_c    = dict.get('SIZE_CRYSTAL')[2]
              n_max  = max(n_a, n_b, n_c)
              print dict.get('NPROD'), dict.get('NCONFIG'), int(dict.get('NPROD')/dict.get('NCONFIG'))
              dict.update({
                   'NMOL'      : prod( dict.get('SIZE_CRYSTAL')),
                   'NATOM_MOL' : sum( 1 for line in open(self._templates_path + dict.get('FILEMOL'))),
                   'NATOMS'    : prod( dict.get('SIZE_CRYSTAL'))* \
                                 sum( 1 for line in open(self._templates_path + dict.get('FILEMOL'))),
                   'NORM_LATTICE': [norm_a, norm_b, norm_c],
                   'LBOXA'       :  10*norm_max*n_max,
                   'LBOXB'       :  10*norm_max*n_max,
                   'LBOXC'       :  10*norm_max*n_max,
                   'RCUT'        :  5*norm_max*n_max,
                   'PRINT'   : int(dict.get('NPROD')/dict.get('NCONFIG'))
                   })

        def write(self):
            self.create_dir()
            self.print_info()
            os.chdir(self._templates_path)
            self._construct_organic_crystal()
            os.chdir(self._bucket_path)

	def _construct_organic_crystal(self):
            molfile = open(self._filemol, 'r')
            crystalfile = open(self._filecrystal, 'w')
            coordfile = open('COORD.temp', 'w')
            lines   = molfile.readlines()
            molfile.close
            for mol0_index in range(self._sizecrystal[0]):
                for mol1_index in range(self._sizecrystal[1]):
                    for mol2_index in range(self._sizecrystal[2]):
                        index3d = [mol0_index+1, mol1_index+1, mol2_index+1] 
                        for line_index in range(len(lines)):
                            l = string.strip(lines[line_index])
                            info = re.split('\s+',l)
                            atom_label = self._choose_atom_label(info[0], i3d = index3d, i3dcharge = self._coordcharge)
                            atom_coord = [float(info[1]), float(info[2]), float(info[3])]
                            vec_shift = mol0_index * array(self._vecta) + \
                                        mol1_index * array(self._vectb) + \
                                        mol2_index * array(self._vectc)
                            result = '%s  %s\n' \
                                 % (atom_label, str(atom_coord + vec_shift).strip('[]') )
                            crystalfile.write(result)
                            coordfile.write(result)
            crystalfile.close
            os.system(' mv %s %s' % (self._filecrystal, self.parcel.path))
            coordfile.close
        
        def create_dir(self):
            os.chdir(self._output_path)
            self.parcel = Dir('CRYSTAL', {})
            self.parcel.mkdir()
            os.chdir(self._bucket_path)


        def _choose_atom_label(self, atom, i3d = [-1, -1, -1], i3dcharge = [-1, -1, -1], imol = -1, icharge = -1):
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


class Dir(object):
      """
      """
      def __init__(self, name, paths = {}):
          self._name = name
          self.path = os.getcwd() + '/' + self._name + '/'
          paths.update({ self._name : self.path})


      def clean(self):
          test = os.path.exists(self.path)
          if (test):
             os.system('rm %s/*.temp ' % self.path)

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


class CP2KFistInput(object):
      """
      """
      def __init__(self, dict, paths, **kwargs):
          self._my_sed_dict = sed_dict
          self._my_sed_dict.update(**kwargs)
          self._my_sed_dict.update(dict)
          self._dict = dict
          self._timestep = dict.get('TIMESTEP')
          self._lengtheq = dict.get('NEQ')
          self._lengthprod = dict.get('NPROD')
          self._sizecrystal = dict.get('SIZE_CRYSTAL')
          self._coordcharge = dict.get('COORD_CHARGE')
          self._mol_name    = dict.get('MOL_NAME')
          self._template_file = dict.get('TEMPLATE_FILE')
          self._filemol = dict.get('FILEMOL')
          self._restraint = 10E-03
          self._norm_lattice = dict.get('NORM_LATTICE')
          self._templates_path = paths.get('templates')
          self._bucket_path = paths.get('bucket')
          self._cp2k_path   = paths.get('cp2k')
 
      def print_info(self):
          print """\n
             TIMESTEP        = %f
             LENGTH          = %d
             PRINT FREQUENCY = %d
             BOX SIZE        = %s
             \n
             """ \
             % ( self._timestep, 
                 self._lengtheq + self._lengthprod,
                 self._dict.get('PRINT'),
                 str(self._sizecrystal).strip('[]') )

      def run(self, ndir):
          self.print_info()
          self.write()
          self.ndir = ndir
          dir = Dir('run-%d' % ndir)
          dir.rm_mkdir()
          os.system('mv %srun.inp %s' % (self._templates_path, dir.path) )
          os.system('cp %s/*psf %s' % ( self._templates_path, dir.path) )
          complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
          print "CP2K STARTS AT: " + complete_time
          dir.chdir()
          os.system( self._cp2k_path + '  run.inp > run.log'   )
          os.chdir(self._bucket_path)
          complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
          print "CP2K FINISHES AT: " + complete_time
          return ndir + 1

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


      def write(self):
          os.chdir(self._templates_path)
          self._psf()
          self._colvar()
          self._constraint()
          input_file = open('run.inp', 'w')
          template_file = open(self._template_file)
          for line in template_file.readlines():
              if 'INCLUDE' in line:
                 result = self._include(line)
              elif 'sed' in line:
                 result = self._sed(line)
              else:
                 result = line
              input_file.write(result)
          input_file.close()
          os.chdir(self._bucket_path)

      def _constraint(self):
          fileout = open('CONSTRAINT.temp', 'w')
          fileout.write('        &CONSTRAINT\n')
          list_mol = [[ i, j, k] for i in range(1, self._sizecrystal[0]+1) 
                                 for j in range(1, self._sizecrystal[1]+1)
                                 for k in range(1, self._sizecrystal[2]+1) ]
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
                         \n     """  \
                         % (  self._restraint,
                              colvar,
                              dot(array(self._norm_lattice), top_vect))
                     fileout.write(result)
          
          fileout.write('        &END CONSTRAINT\n')
          fileout.close()

      def _colvar(self):
          fileout = open('COLVAR.temp', 'w')
          list_mol = [[ i, j, k] for i in range(1, self._sizecrystal[0]+1) 
                                 for j in range(1, self._sizecrystal[1]+1)
                                 for k in range(1, self._sizecrystal[2]+1) ]
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
                """  %( \
                      self._list_carbon(mol),
                      self._list_carbon(mol2))
                     fileout.write(result)
          fileout.close()
                   
      def _list_carbon(self, mol):
            molfile = open(self._filemol, 'r')
            lines = molfile.readlines()
            list_carbon = []
            num_mol    = (mol[0] - 1 )*self._sizecrystal[1]*self._sizecrystal[2] + \
                         (mol[1] - 1 )*self._sizecrystal[2] + \
                          mol[2] - 1 
            for i, line in enumerate(lines):
                if ('C' in line):
                   list_carbon.append(i+1 + (num_mol)*(len(lines)))
            return '\t'.join(map(str, list_carbon))      
 

      def _include(self, line):
          word_list = line.split()
          if word_list[0] != 'INCLUDE':
             print 'PROBLEM IN INCLUDE FUNCTION'
             sys.exit()
          file_to_include = open(word_list[1])
          result = ''
          for my_line in file_to_include.readlines():
              result = result + my_line 
          return result 
          

      def _psf(self):
          fileout = open('PSF.temp', 'w')
          number_mol = prod(self._sizecrystal)
          pos_mol    = (float(self._coordcharge[0]) - 1 )*float(self._sizecrystal[1]*self._sizecrystal[2]) + \
                       (float(self._coordcharge[1]) - 1 )*float(self._sizecrystal[2]) + \
                        float(self._coordcharge[2]) 
          for mol in range(int(number_mol)):
             if ( (mol+1) == pos_mol):
                name = self._mol_name + '_CHARGE.psf'
             else:
                name = self._mol_name + '_NEUTRE.psf'
             result="""\n                               
                    &MOLECULE
                        NMOL              1
                        CONN_FILE_NAME    ./%s
                        CONN_FILE_FORMAT  PSF
                   &END MOLECULE\n """ % name
             fileout.write(result)
          fileout.close()


class FSSHParcel(object):
      """
      """
      def __init__(self, dict, paths):
          self._nmol = dict.get('NMOL')
          self._natom_mol = dict.get('NATOM_MOL')
          self._timestep = dict.get('TIMESTEP')
          self._lengtheq = dict.get('NEQ')
          self._nprod = dict.get('NPROD')
          self._printfrq = dict.get('PRINT')
          self._sizecrystal = dict.get('SIZE_CRYSTAL')
          self._coordcharge = dict.get('COORD_CHARGE')
          self._mol_name    = dict.get('MOL_NAME')
          self._template_file = dict.get('TEMPLATE_FILE')
          self._filemol = dict.get('FILEMOL')
          self._restraint = 10E-03
          self._templates_path = paths.get('templates')
          self._bucket_path = paths.get('bucket')
          self._output_path = paths.get('output')
          self._generator_path   = paths.get('generator')

      def gather(self, ndir):
          self.create()
          self.gather_vel_coord(ndir)
          self.gather_templates_generator()    
          self.prepare_input()

      def prepare_input(self):
          os.chdir(self.parcel.path)
          toname = 'TONAME'
          my_input = open('%s.input' % toname, 'w' )
          result = """
KIND_RUN  FSSH_OS
TEMPLATE_FILE FSSH_CORE.template
FORCEFIELD_FILE FSSH_FF.template
NUMBER_INIT %d
MOL_NAME %s
NATOMS   %d
NATOM_MOL   %d
VECTA %s
VECTB %s
VECTC %s
SIZE_CRYSTAL %s
COORD_CHARGE %s
STEPS        1
PRINTFRQ     1
TEST         NO
        """ % (
          sed_dict.get('NCONFIG', '!!!'),
          sed_dict.get('MOL_NAME', '!!!'),
          sed_dict.get('NATOMS', '!!!'),
          sed_dict.get('NATOM_MOL', '!!!'),
          '   '.join(map(str, sed_dict.get('VECTA', '!!!'))),
          '   '.join(map(str, sed_dict.get('VECTB', '!!!'))),
          '   '.join(map(str, sed_dict.get('VECTC', '!!!'))),
          '   '.join(map(str, sed_dict.get('SIZE_CRYSTAL', '!!!'))),
          '   '.join(map(str, sed_dict.get('COORD_CHARGE', '!!!'))))
          my_input.write(result)
          my_input.close()
         
              

      def gather_vel_coord(self, ndir):
          os.chdir('run-%d' % ndir)
          pos_mol    = (self._coordcharge[0] - 1 )*(self._sizecrystal[1]*self._sizecrystal[2]) + \
                       (self._coordcharge[1] - 1 )*(self._sizecrystal[2]) + \
                       (self._coordcharge[2]) - 1
          for iprop in ['vel','pos']:
              filein = open('run-' + iprop +'-1.xyz', 'r')
              lines = filein.readlines()
              iconfig = 0
              for istep in range(0, self._nprod, self._printfrq):
                  iconfig = iconfig + 1
                  filename = iprop + '-' + str(iconfig) + '.init'
                  fileout = open( filename, 'w')
                  for imol in range(self._nmol):
                         for iatom in range(self._natom_mol):
                             index = (istep/self._printfrq)*(self._nmol*self._natom_mol+2) \
                                   + 2                    \
                                   + imol*self._natom_mol           \
                                   + iatom
                             l = string.strip(lines[index])
                             info = re.split('\s+',l)
                             atom_label = self._choose_atom_label(info[0], imol = imol, icharge = pos_mol)
                             atom_xyz = [float(info[1]), float(info[2]), float(info[3])]
                             result = '%s  %s\n' \
                                 % (atom_label, str(atom_xyz).strip('[]') )
                             fileout.write(result)
                  fileout.close()
                  os.system(' mv %s %s' % ( filename, self.initial.path) )
          os.chdir(self._bucket_path)                

      def gather_templates_generator(self):
          os.system('cp %s/*  %s' % (self._templates_path, self.templates.path))
          os.system('cp %s/*  %s' % (self._generator_path, self.generator.path))

      def create(self):
          os.chdir(self._output_path)
          self.parcel = Dir('TONAME')
          self.parcel.mkdir()
          self.parcel.chdir()
          self.generator = Dir('generator')
          self.generator.mkdir()
          self.templates = Dir('templates')
          self.templates.mkdir()
          self.initial = Dir('initial')
          self.initial.mkdir()
          os.chdir(self._bucket_path)

      def _choose_atom_label(self, atom, i3d = [-1, -1, -1], i3dcharge = [-1, -1, -1], imol = -1, icharge = -1):
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


class CP2KFSSH(object):
      """
      """
      def __init__(self, dict, paths, **kwargs):
          self._my_sed_dict = sed_dict
          self._my_sed_dict.update(**kwargs)
          self._my_sed_dict.update(dict)
          self._init         = self._my_sed_dict.get('INIT')
          self._timestep = dict.get('TIMESTEP')
          self._printfrq = dict.get('PRINTFRQ')
          self._sizecrystal = dict.get('SIZE_CRYSTAL')
          self._coordcharge = dict.get('COORD_CHARGE')
          self._mol_name    = dict.get('MOL_NAME')
          self._template_file = dict.get('TEMPLATE_FILE')
          self._forcefield_file = dict.get('FORCEFIELD_FILE')
          self._filemol = 'COORD.temp'
          self._restraint = 10E-03
          self._norm_lattice = dict.get('NORM_LATTICE')
          self._templates_path = paths.get('templates')
          self._bucket_path = paths.get('bucket')
          self._cp2k_path   = paths.get('cp2k')
          self._initial_path   = paths.get('initial')
          self._natom_mol      = dict.get('NATOM_MOL')
 
      def print_info(self):
           print "Hey"
#          print """\n
#             TIMESTEP        = %f
#             LENGTH          = %d
#             PRINT FREQUENCY = %d
#             BOX SIZE        = %s
#             \n
#             """ \
#             % ( self._timestep, 
#                 self._lengtheq + self._lengthprod,
#                 self._printfrq,
#                 str(self._sizecrystal).strip('[]') )

      def run(self, ndir):
          self.ndir = ndir
          self.print_info()
          self.write()
          dir = Dir('run-%d' % ndir)
          dir.rm_mkdir()
          os.system('mv %srun.inp %s' % (self.tmp.path, dir.path) )
          os.system('cp %s/*psf %s' % ( self.tmp.path, dir.path) )
          complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
          print "CP2K STARTS AT: " + complete_time
          dir.chdir()
          val = os.system( self._cp2k_path + '  run.inp > run.log'   )
          os.chdir(self._bucket_path)
          complete_time = time.strftime("%y%m%d%H%M%S", time.localtime())
          print "CP2K FINISHES AT: " + complete_time
          if val == 0:
             return ndir + 1
          else:
             return -1

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

      def create_tmp(self):
          self.tmp = Dir('tmp')
          self.tmp.rm_mkdir()
          os.system('cp %s/*.psf %s' % (self._templates_path, self.tmp.path) )
          os.system('cp %s/*.inc %s' % (self._templates_path, self.tmp.path) )
          os.system('cp %s/FSSH* %s' % (self._templates_path, self.tmp.path) )
          os.system('cp %s/pos-%d.init %s/COORD.temp' % (self._initial_path, self._init, self.tmp.path) )
          os.system('cp %s/vel-%d.init %s/VELOC.temp' % (self._initial_path, self._init, self.tmp.path) )

      def _complete_dict(self):
          dict = self._my_sed_dict
          if (dict) is None:
              print "Error in input.complete"
              sys.exit()
          else:
              norm_a = norm(array(dict.get('VECTA')))
              norm_b = norm(array(dict.get('VECTB')))
              norm_c = norm(array(dict.get('VECTC')))
              n_a    = dict.get('SIZE_CRYSTAL')[0]
              n_b    = dict.get('SIZE_CRYSTAL')[1]
              n_c    = dict.get('SIZE_CRYSTAL')[2]
              n_max  = max(n_a, n_b, n_c)
              norm_max = max(norm_a, norm_b, norm_c)
              dict.update({
                   'NMOL'      : prod( dict.get('SIZE_CRYSTAL')),
                   'NDIABAT'     :  prod(dict.get('SIZE_CRYSTAL'))*dict.get('NORBITALS'),
                   'NORM_LATTICE': [norm_a, norm_b, norm_c],
                   'LBOXA'       :  10*norm_max*n_max,
                   'LBOXB'       :  10*norm_max*n_max,
                   'LBOXC'       :  10*norm_max*n_max,
                   'RCUT'        :  5*norm_max*n_max
                   })
              dict.update({
                   'FORCE_EVAL_ORDER' : '  '.join(map(str, range(1,dict.get('NDIABAT')+2)) )
                         })
          self._nmol = dict.get('NMOL')
          self._norm_lattice = dict.get('NORM_LATTICE')

      def write(self):
          self.create_tmp()
          self._complete_dict()
          os.chdir(self.tmp.path)
          for mol in range(1, self._nmol+1):
              self._psf(mol)
          self._colvar()
          self._constraint()
          self._aom()
          force_field = open('FORCEFIELD.tmp', 'w')
          pre_ff = open(self._forcefield_file)
          result = self._amend(pre_ff, self._nmol)
          force_field.write(result)
          force_field.close()

          input_file = open('run.inp', 'w')
          template_file = open(self._template_file)
          result = self._amend(template_file)
          input_file.write(result)
          input_file.close()
          os.chdir(self._bucket_path)

      def _aom(self):
          for mol in range(self._nmol):
              os.system('cat %s_AOM.inc >> AOM_COEFF.temp' % self._mol_name )

      def _amend(self, file, number = 1):
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

      def _constraint(self):
          fileout = open('CONSTRAINT.temp', 'w')
          fileout.write('        &CONSTRAINT\n')
          list_mol = [[ i, j, k] for i in range(1, self._sizecrystal[0]+1) 
                                 for j in range(1, self._sizecrystal[1]+1)
                                 for k in range(1, self._sizecrystal[2]+1) ]
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
                         \n     """  \
                         % (  self._restraint,
                              colvar,
                              dot(array(self._norm_lattice), top_vect))
                     fileout.write(result)
          
          fileout.write('        &END CONSTRAINT\n')
          fileout.close()

      def _colvar(self):
          fileout = open('COLVAR.temp', 'w')
          list_mol = [[ i, j, k] for i in range(1, self._sizecrystal[0]+1) 
                                 for j in range(1, self._sizecrystal[1]+1)
                                 for k in range(1, self._sizecrystal[2]+1) ]
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
                """  %( \
                      self._list_carbon(mol),
                      self._list_carbon(mol2))
                     fileout.write(result)
          fileout.close()
                   
      def _list_carbon(self, mol):
            molfile = open(self._filemol, 'r')
            lines = molfile.readlines()
            list_carbon = []
            num_mol    = (mol[0] - 1 )*self._sizecrystal[1]*self._sizecrystal[2] + \
                         (mol[1] - 1 )*self._sizecrystal[2] + \
                          mol[2] - 1 
            for i in range(self._natom_mol):
                line = lines[i]
                if ('C' in line):
                   list_carbon.append(i+1 + (num_mol)*self._natom_mol)
            return '\t'.join(map(str, list_carbon))      
 

      def _include(self, line, mol = 0):
          word_list = line.split()
          if word_list[0] != 'INCLUDE':
             print 'PROBLEM IN INCLUDE FUNCTION'
             sys.exit()
          if word_list[1] == 'SPECIAL':
             file_to_include = open(word_list[2] + '-' + str(mol) + '.temp')
          else:
             file_to_include = open(word_list[1])
          result = ''
          for my_line in file_to_include.readlines():
              result = result + my_line 
          return result 
          

      def _psf(self, pos_mol):
          fileout = open('PSF-%d.temp' % pos_mol, 'w')
          number_mol = prod(self._sizecrystal)
          #pos_mol    = (float(self._coordcharge[0]) - 1 )*float(self._sizecrystal[1]*self._sizecrystal[2]) + \
          #             (float(self._coordcharge[1]) - 1 )*float(self._sizecrystal[2]) + \
          #              float(self._coordcharge[2]) 
          for mol in range(int(number_mol)):
             if ( (mol+1) == pos_mol):
                name = self._mol_name + '_CHARGE.psf'
             else:
                name = self._mol_name + '_NEUTRE.psf'
             result="""\n                               
                    &MOLECULE
                        NMOL              1
                        CONN_FILE_NAME    ./%s
                        CONN_FILE_FORMAT  PSF
                   &END MOLECULE\n """ % name
             fileout.write(result)
          fileout.close()



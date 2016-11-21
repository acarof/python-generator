#!/usr/bin/python

import string, re, struct, sys, math, os, time
import gentool
import numpy

try:
         finput_name = sys.argv[1]
except: 
          raise SystemExit

# SET_UP THE DIRECTORY, CHECK ANY SUBDIR IS PRESENT
gentool.mkdir('CRYSTAL')
gentool.checkdir('MOLECULE')
gentool.checkdir('TEMPLATE_CP2K')


# DEFINE PATH/EXE VARIABLES
pcentral = os.getcwd()
ptemplate = pcentral + '/TEMPLATE_CP2K/'
pexe = '/scratch/grudorff/antoine/bin'
xcp2k = pexe + '/cp2k.sopt'

# OPEN AND READ THE INPUT FILE
finput = open( finput_name, 'r')
lines  = finput.readlines()
finput.close()

fmol_name = lines[0].rstrip('\r\n')
l           = string.strip(lines[1])
info        = re.split('\s+', l)
veca_lattice = [float(info[0]), float(info[1]), float(info[2])]
l           = string.strip(lines[2])
info        = re.split('\s+', l)
vecb_lattice = [float(info[0]), float(info[1]), float(info[2])]
l           = string.strip(lines[3])
info        = re.split('\s+', l)
vecc_lattice = [float(info[0]), float(info[1]), float(info[2])]
l           = string.strip(lines[4])
info        = re.split('\s+', l)
size_crystal = [int(info[0]), int(info[1]), int(info[2])]
l             = string.strip(lines[5])
info          = re.split('\s+', l)
icoord_charge = [int(info[0]), int(info[1]), int(info[2])]
fcrystal_name = lines[6].rstrip('\r\n')
molecule      = lines[7].rstrip('\r\n')
timestep      = float(lines[8].rstrip('\r\n'))
nsteps_eq     = int(lines[9].rstrip('\r\n'))
nsteps_prod   = int(lines[10].rstrip('\r\n'))
nconfig       = int(lines[11].rstrip('\r\n'))
l             = string.strip(lines[12])
info          = re.split('\s+', l)
sizebox       = [float(info[0]), float(info[1]), float(info[2])]
test     = bool(lines[13].rstrip('\r\n'))


# DEFINE SOME VARIABLE
is_a_test = test == 'TEST'
pre_bucket = 'DIABAT_OS_' +  fcrystal_name
bucket = gentool.name_bucket(pre_bucket, is_a_test)
print_frq = int(math.floor(nsteps_prod / nconfig))
nmol      = numpy.prod(size_crystal)
natom_mol = sum( 1 for line in open('MOLECULE/' + fmol_name))
natoms    = nmol*natom_mol
nsteps    = nsteps_eq + nsteps_prod


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print """!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
START THE DIABAT OF SAMPLING OF AN ORGANIC CRYSTAL OF %s
BUCKET = %s
""" % ( molecule, bucket)

print """!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1. CONSTRUCT THE ORGANIC CRYSTAL.
MOLECULE FILE   = %s
CRYSTAL LATTICE a = %f, %f, %f
CRYSTAL LATTICE b = %f, %f, %f
CRYSTAL LATTICE c = %f, %f, %f
CRYSTAL SIZE    = %d, %d, %d
COORD CHARGE    = %d, %d, %d
CRYSTAL FILE    = %s.coord
""" \
% (fmol_name, veca_lattice[0], veca_lattice[1], veca_lattice[2], \
vecb_lattice[0], vecb_lattice[1], vecb_lattice[2], \
vecc_lattice[0], vecc_lattice[1], vecc_lattice[2], \
size_crystal[0], size_crystal[1], size_crystal[2], \
   icoord_charge[0], icoord_charge[1], icoord_charge[2], fcrystal_name)

os.chdir('MOLECULE')
gentool.construct_organic_crystal( fmol_name, veca_lattice, vecb_lattice, vecc_lattice, size_crystal, \
         icoord_charge, fcrystal_name + '.coord')
if not is_a_test:
   os.chdir('../CRYSTAL')
   gentool.mkdir_bucket(bucket)
   os.system('cp ../MOLECULE/%s.coord  %s/' % (fcrystal_name, bucket) )
   os.system('mv ../MOLECULE/%s.coord  ..' %  (fcrystal_name) )
else:
   os.system('mv %s.coord  ..' %  (fcrystal_name) )
   
os.chdir('..')



print """!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
2. RUN CP2K
TIMESTEP        = %d
LENGTH          = %d
PRINT FREQUENCY = %d
BOX SIZE        = %f, %f, %f
""" % (timestep, nsteps_eq + nsteps_prod, print_frq, sizebox[0], sizebox[1], sizebox[2])

number_dir = -1
os.chdir('..')
gentool.mkdir_bucket(bucket)
os.chdir(bucket)
gentool.mkdir('run0')
os.chdir('run0')
os.system('cp %s/FIST_TEMPLATE .' % ptemplate)
os.system('cp %s/%s* .'           % (ptemplate, molecule))
os.system('mv %s/%s.coord COORD.init'   % (pcentral, fcrystal_name))
gentool.input_fist(molecule, timestep, nsteps, print_frq, natoms, sizebox)
gentool.psf( size_crystal, icoord_charge, molecule)
gentool.constraint( nmol, natom_mol, molecule)
os.system(xcp2k + ' run.inp > run.out 2>   run.err')



print """
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
3. GATHER VELOCITIES AND COORDINATE
"""


gentool.gather_config(molecule, nmol, natom_mol, icoord_charge, size_crystal, nsteps_eq, nsteps_prod, print_frq)
os.system('cp vel* pos* %s/CRYSTAL/%s ' %(pcentral, bucket) )
os.chdir(pcentral)


print """
SAMPLING OF THE DIABATIC OS IS OVER!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""

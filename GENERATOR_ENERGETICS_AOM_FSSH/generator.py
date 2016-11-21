#!/usr/bin/python

import string, re, struct, sys, math, os, time
import gentool

#try : 
#         DATE         = sys.argv[1]
#         CP2K_VERSION = sys.argv[2]
#except :
#         raise SystemExit

#MODULEPATH = '/scratch/grudorff/modulefiles module load mpich'
#PWD = os.getcwd()
#PATH_EXE = '/scratch/grudorff/antoine/bin'
#PATH_TEMPLATE = PWD  + '/TEMPLATE/'
#PATH_CONFIG = PWD + '/CONFIG_INITIAL'
#exe_cp2k = PATH_EXE + '/cp2k_nonadiabatic_scutum_160427.sopt'

#BUCKET='ENERGIES_TEST_%s_%s' % (DATE, CP2K_VERSION)
#os.chdir('..')
#print os.getcwd()
#test_bucket = os.path.exists(BUCKET)
#if (test_bucket):
#   print "THE DIRECTORY EXISTS ALREADY!"
#   sys.exit()
#else:
#   os.mkdir(BUCKET)
#   os.chdir(BUCKET)


bucket = 'ENERGIES'
gentool.create_bucket(bucket)
number_dir = -1


print """!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
STUDY BORN_OPPENHEIMER FOR DIFFERENT SCALING FACTOR (COUPLING)
"""
scalingList = [0.0065, 0.0130, 0.0325, 0.0652, 0.1956, 0.3912, 0.9780]  # in eV
diabatList  = [1, 2]
methodList  = ['BORN_OPPENHEIMER']
name_coord = 'COORD'
name_veloc = 'VELOC.init'
configList = [1]
stepsList    = [ 10 ]
configList = [1]
number_dir  = gentool.create_and_run(SCALINGLIST = scalingList, DIABATLIST = diabatList, METHODLIST = methodList, \
                             number_dir = number_dir, STEPSLIST = stepsList, NAME_COORD = name_coord,             \
                             NAME_VELOC = name_veloc, CONFIGLIST = configList)



print """!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
STUDY SURFACE HOPPING FOR DIFFERENT SCALING FACTOR
"""
scalingList = [0.0065, 0.0130, 0.0325, 0.0652, 0.1956, 0.3912, 0.9780]  # in eV
methodList  = ['FSSH']
name_coord = 'COORD'
name_veloc = 'VELOC.init'
configList = [1]
stepsList    = [ 10 ]
number_dir  = gentool.create_and_run(SCALINGLIST = scalingList,  METHODLIST = methodList, \
                             number_dir = number_dir, STEPSLIST = stepsList, NAME_COORD = name_coord,             \
                             NAME_VELOC = name_veloc, CONFIGLIST = configList)


print """
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
STUDY SURFACE HOPPING FOR DIFFERENT TIMESTEP 
"""
timestepList = [0.05, 0.1, 0.5, 1]
stepsList    = [ 10 , 10 , 10 ,  10]
scalingList =  [0.0065]
methodList  = ['FSSH']
name_coord = 'COORD'
name_veloc = 'VELOC.init'
configList = [1]
number_dir  = gentool.create_and_run(SCALINGLIST = scalingList, METHODLIST = methodList, \
                             number_dir = number_dir, STEPSLIST = stepsList, NAME_COORD = name_coord,             \
                             NAME_VELOC = name_veloc, CONFIGLIST = configList,   \
                             TIMESTEPLIST = timestepList)









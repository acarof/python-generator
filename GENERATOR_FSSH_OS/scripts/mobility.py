# standard modules
import string, re, struct, sys, math, os, time
import numpy as np
import importlib, imp
from multiprocessing import Pool, cpu_count
from functools import partial
from datetime import datetime

# custom modudules
from utils_scripts import *
#from marcus import *
from utils_analyse import *

scripts = 'mobility'
keywords = ['TEMPERATURE']
dict_properties = {
    'Block-runs-average' : ['MSD']
}
number_blocks = 2

#FIND THE TITLE
if 'GENERATOR_FSSH_OS' in os.getcwd():
    title = 'TEST'
else:
    name_bucket = os.getcwd().split('/')[-1]
    short_time = time.strftime("%y%m%d-%H%M", time.localtime())
    title = '%s-%s' % (name_bucket, short_time,)
dataname = 'data-%s-%s' % (scripts, title)
if not os.path.isdir(dataname):
    os.mkdir(dataname)


# CREATE LIST OF RUNS
print "Start to build run_dir at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")
dirlist = os.listdir('.')
run_dict = {}
os.system('cd ..')
for i, directory in enumerate(dirlist):
    if 'run-' in directory and 'per' not in directory:
        dir = FSSHRun(directory)
        keys = dir.get_input_key(keywords)
        if run_dict.get(keys ) is None:
            run_dict[keys] = []
        run_dict[keys].append(directory)
filename = print_run_dict(run_dict, keywords)
os.system('mv %s %s' % (filename, dataname))
print "Finish to build run_dir at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")


# PARALLEL OR SERIAL CALCULATION
def super_analyse(tuple):
    return analyse_properties(tuple, run_dict, dict_properties, number_blocks=number_blocks)
try:
    nworker = int( sys.argv[1] )
except:
    nworker = cpu_count()
if nworker == 0:
    results = []
    for keys in run_dict.keys():
        results.append( super_analyse(keys) )
else:
    pool = Pool(nworker)
    results = pool.map(super_analyse, run_dict.keys())


# PRINT THE RESULTS
results_dict = {}
filetuple = open('List-tuple.dat', 'w')
line = '%s\n' % ('  '.join(keywords))
filetuple.write(line)
for tuple, result in zip( run_dict.keys(), results):
    filetuple.write('%s\n' % '  '.join(tuple))
    results_dict[tuple] = result
    properties = results_dict[tuple]
    for property in dict_properties.get('Block-runs-average', []):
        new_list = []
        for element in properties[property]:
            new_list.append(average_dict(element, properties['Length-block']))
        filename = print_list_dict(new_list, property, tuple)
        os.system('mv %s %s' % (filename, dataname))
filetuple.close()
os.system('mv List-tuple.dat %s' % dataname )
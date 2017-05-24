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

scripts = 'extract-scaling-deco'
keywords = ['SCALING_FACTOR']

dict_properties = {
    #  'Detailled' : [ 'Adiabatic-populations','Surface-populations','Delta_E'],
    'Runs-average' : [],
    'Block-runs-average' : ['Adiabatic-populations','Surface-populations','Delta_E'],
    #    'Specific' : ['FSSH', 'Detailed-FSSH'],
    'Specific' : ['FSSH', 'Detailed-FSSH'],
    #'Initial' : ['Delta_E'],
    'Initial' : [],
    #'Mean' : ['Temperature','Couplings', 'Total-energy']
    'Mean' : ['Temperature','Couplings', 'Total-energy']
}

#dict_properties = {
#    'Initial' : ['Couplings']
#}

#  INPUT INFO
try:
    info = sys.argv[1]
except:
    info = None
previous_name = None
unique_run = None
if info is not None:
    if 'run-' in info:
        unique_run = info.split('/')[0]
    elif 'data-' in info:
        previous_name = info



#FIND THE TITLE
if 'GENERATOR_GLOBAL' in os.getcwd():
    dirlist = os.listdir('.')
    title = 'TEST'
else:
    name_bucket = os.getcwd().split('/')[-1]
    short_time = time.strftime("%y%m%d%H%M", time.localtime())
    title = '%s-%s' % (name_bucket, short_time,)



#FIND THE NAME OF THE DATA DIRECTORY
if unique_run is None:
    dirlist = os.listdir('.')
    if previous_name is None:
        dataname = 'data-%s-%s' % (scripts, title)
        if not os.path.isdir(dataname):
            os.mkdir(dataname)
    else:
        dataname = previous_name
else:
    dirlist = [unique_run]
    dataname = 'one-%s-%s-%s' % (unique_run, scripts, title)
    if not os.path.isdir(dataname):
        os.mkdir(dataname)



# CREATE LIST OF RUNS
print "Start to build run_dir at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")
run_dict = {}
os.system('cd ..')
if previous_name is None:
    for i, directory in enumerate(dirlist):
        if 'run-' in directory and 'per' not in directory:
            dir = FSSHRun(directory)
            keys = dir.get_input_key(keywords)
            if run_dict.get(keys ) is None:
                run_dict[keys] = []
            run_dict[keys].append(directory)
    filename = print_run_dict(run_dict, keywords)
    os.system('mv %s %s' % (filename, dataname))
else:
    file_listdir = open('%s/List-run.dat' % dataname)
    for line in file_listdir.readlines():
        list = re.split(r'(\(.*\))(.*)', line)
        tuple_ = tuple(re.findall(r'([^,|^)|^(|^ ]+)', list[1].replace("'", "")))
        runs = list[2].split()
        run_dict[tuple_] = runs
print "Finish to build run_dir at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")



# PARALLEL OR SERIAL CALCULATION
def super_analyse(tuple):
    return analyse_properties(tuple, run_dict, dict_properties)
try:
    nworker = int( sys.argv[1] )
except:
    nworker = -1
if nworker == -1:
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
    for property in dict_properties.get('Mean', []):
        filename = create_file(property, title, properties[property + 'info'], 'Mean', tuple = tuple)
        os.system('mv %s %s' % (filename, dataname))
    for property in dict_properties.get('Runs-average', []):
        properties[property] = average_dict(properties[property], properties['Number runs'])
        filename = print_dict( properties[property], property, tuple  )
        os.system('mv %s %s' % (filename, dataname))
    for property in dict_properties['Block-runs-average']:
        new_list = []
        for element in properties[property]:
            new_list.append(average_dict(element, properties['Length-block']))
        filename = print_list_dict(new_list, property, tuple)
        os.system('mv %s %s' % (filename, dataname))
    for property in dict_properties.get('Specific', []):
        filename = create_file(property, title, properties[property + 'info'], 'Spec', tuple = tuple)
        os.system('mv %s %s' % (filename, dataname))
    for property in dict_properties.get('Initial', []):
        filename = create_file(property, title, properties[property + 'info'], 'Initial', tuple = tuple)
        os.system('mv %s %s' % (filename, dataname))
filetuple.close()
os.system('mv List-tuple.dat %s' % dataname )



print "End of the analysis at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")



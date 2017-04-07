# standard modules
import string, re, struct, sys, math, os, time
import numpy as np
import importlib, imp
#import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
#from operator import itemgetter
from multiprocessing import Pool, cpu_count
#from scipy.integrate import quad

# custom modudules
from utils_scripts import *
from marcus import *
from datetime import datetime

scripts = 'extract-scaling-deco'
keywords = ['SCALING_FACTOR','DECOHERENCE_CORRECTIONS']

detail_properties = ['Surface-populations', 'Adiabatic-populations']
#ratio_properties = ['Internal consistency ratio']
#histo_properties = ['Delta_E', 'Couplings', 'Populations']
histo_properties = ['Delta_E']
mean_properties = ['Temperature','Couplings']
specific_properties = ['FSSH', 'Detailed-FSSH']
total_properties = detail_properties + mean_properties + specific_properties + histo_properties

if 'GENERATOR_GLOBAL' in os.getcwd():
    dirlist = os.listdir('.')
    title = 'TEST'
else:
    dirlist = os.listdir('.')
    name_bucket = os.getcwd().split('/')[-1]
    short_time = time.strftime("%y%m%d%H%M", time.localtime())
    title = '%s-%s' % (name_bucket, short_time,)
dataname = 'data-%s-%s' % (scripts, title)
if not os.path.isdir(dataname):
    os.mkdir(dataname)

bin = 2
bin_histo = 50

reorga = 0.300
free_energy = 0.00
nadiab = 2


def sum_two_dict( dict1, dict2):
    result = {}
    if dict1 is None:
        result = dict2
    else:
        for key in dict1:
            value = np.array(dict1[key]) + np.array(dict2[key])
            result[key] = value
    return result


def average_dict(dict1, number):
    result = {}
    for key in dict1:
        result[key] = np.array( dict1[key] ) / number
    return result


def print_dict( dict_, property, tuple):
    filename = property + '-' + '-'.join(tuple) + '.dat'
    file = open(filename, 'w')
    for time in sorted(dict_):
        line = '%f  %s\n' % (time, '   '.join(map(str, dict_[time])) )
        file.write(line)
    file.close()
    return filename


def analyse_properties(tuple):
    real_start_time = datetime.now()
    start_time = datetime.strftime( real_start_time, "%Y %m %d %H:%M:%S ")
    print "One worker for: %s starts at %s" % (tuple, start_time)
    list_dir = run_dict[tuple]
    properties_dict = {}
    for directory in list_dir:
        os.chdir(directory)
        dir = FSSHRun(directory)
        for property in (total_properties):
            prop = dir.extract(property)
            if property in detail_properties:
                properties_dict[property] = sum_two_dict(properties_dict.get(property), prop)
            elif property in mean_properties:
                list = statistics( prop )
                if properties_dict.get(property) is None:
                    properties_dict[property] = []
                    properties_dict[property + 'info'] = ''
                properties_dict[property].append(list)
                line = '%s      %s\n' % ( directory, '    '.join(map(str, list )))
                properties_dict[property + 'info'] += line
            elif property in specific_properties:
                list = prop
                if properties_dict.get(property) is None:
                    properties_dict[property] = []
                    properties_dict[property + 'info'] = ''
                line = '%s      %s\n' % (directory, '    '.join(map(str, list)))
                properties_dict[property + 'info'] += line
        os.chdir('..')
    properties_dict['Number runs'] = len(list_dir)
    real_end_time = datetime.now()
    end_time = datetime.strftime( real_end_time, "%Y %m %d %H:%M:%S ")
    print "One worker for: %s ends at %s" % (tuple, end_time)
    time_diff = real_end_time - real_start_time
    hours, remainder = divmod(time_diff.total_seconds(), 3600)
    minutes, seconds = divmod(remainder, 60)
    print "The worker %s lasted: %s hours %s minutes %s seconds" % (tuple, hours, minutes, seconds)
    return properties_dict


def print_run_dict(run_dict, keywords):
    filename = 'List-run.dat'
    file = open(filename, 'w')
    for tuple in run_dict:
        line = '%s %s\n' % ( tuple, '  '.join(run_dict[tuple]) )
        file.write(line)
    file.close()
    os.system('mv %s %s' % (filename, dataname))


# CREATE LIST OF RUNS
print "Start to build run_dir at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")
run_dict = {}
os.system('cd ..')
for i, directory in enumerate(dirlist):
    if 'run-' in directory and 'per' not in directory:
        os.chdir(directory)
        dir = FSSHRun(directory)
        keys = dir.get_input_key(keywords)
        if run_dict.get(keys ) is None:
            run_dict[keys] = []
        run_dict[keys].append(directory)
        os.chdir('..')
print_run_dict(run_dict, keywords)
print "Finish to build run_dir at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")


# PARALLEL OR SERIAL CALCULATION
try:
    nworker = int( sys.argv[1] )
except:
    nworker = -1
if nworker == -1:
    nworker = cpu_count()
if nworker == 0:
    results = []
    for keys in run_dict.keys():
        results.append( analyse_properties(keys) )
else:
    pool = Pool(nworker)
    results = pool.map(analyse_properties, run_dict.keys())


# PRINT THE RESULTS
results_dict = {}
filetuple = open('List-tuple.dat', 'w')
line = '%s\n' % ('  '.join(keywords))
filetuple.write(line)
for tuple, result in zip( run_dict.keys(), results):
    filetuple.write('%s\n' % '  '.join(tuple))
    results_dict[tuple] = result
    properties = results_dict[tuple]
    for property in mean_properties :
        filename = create_file(property, title, properties[property + 'info'], 'Mean', tuple = tuple)
        os.system('mv %s %s' % (filename, dataname))
    for property in detail_properties:
        properties[property] = average_dict(properties[property], properties['Number runs'])
        filename = print_dict( properties[property], property, tuple  )
        os.system('mv %s %s' % (filename, dataname))
    for property in specific_properties:
        filename = create_file(property, title, properties[property + 'info'], 'Spec', tuple = tuple)
        os.system('mv %s %s' % (filename, dataname))
filetuple.close()
os.system('mv List-tuple.dat %s' % dataname )



print "End of the analysis at %s" %\
    datetime.strftime(datetime.now(), "%Y %m %d %H:%M:%S ")



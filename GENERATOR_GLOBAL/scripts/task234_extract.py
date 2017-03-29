import string, re, struct, sys, math, os, time
import numpy as np
import importlib, imp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from operator import itemgetter
from multiprocessing import Pool
from scipy.integrate import quad

from utils_scripts import *
from marcus import *
from datetime import datetime

scripts = 'task234_parallel'

detail_properties = ['Surface-populations', 'Adiabatic-populations']
#ratio_properties = ['Internal consistency ratio']
#histo_properties = ['Delta_E', 'Couplings', 'Populations']
histo_properties = ['Delta_E']
mean_properties = ['Temperature','Couplings']
#specific_properties = ['FSSH']
specific_properties = []
total_properties = detail_properties + mean_properties + specific_properties + histo_properties
test = True

if test:
    #dirlist = ['run-%d' % i for i in range(10)]
    dirlist = os.listdir('.')
    title = 'TEST'
else:
    dirlist = os.listdir('.')
    name_bucket = os.getcwd().split('/')[-1]
    short_time = time.strftime("%y%m%d%H%M", time.localtime())
    title = '%s-%s' % (name_bucket, short_time,)


bin = 2
bin_histo = 50
reorga = 0.300
free_energy = 0.00
nadiab = 2


dataname = 'data-%s-%s' % (scripts, title)
if not os.path.isdir(dataname):
    os.mkdir(dataname)



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


this_time = datetime.now()
this_time_str = datetime.strftime(this_time, "%Y %m %d %H:%M:%S ")
print "Start to build run_dir at %s" % (this_time_str)

run_dict = {}
os.system('cd ..')
for i, directory in enumerate(dirlist):
    if 'run-' in directory and 'per' not in directory:
        os.chdir(directory)
        dir = FSSHRun(directory)
        ( scaling, reversal ) = dir.get_input_key(['SCALING_FACTOR','METHOD_REVERSAL'])
        if run_dict.get( ( scaling, reversal ) ) is None:
            run_dict[( scaling, reversal )] = []
        run_dict[( scaling, reversal )].append(directory)
        os.chdir('..')


def print_run_dict(run_dict):
    filename = 'List-run.dat'
    file = open(filename, 'w')
    for tuple in run_dict:
        line = '%s %s\n' % ( tuple, '  '.join(run_dict[tuple]) )
        file.write(line)
    file.close()
    os.system('mv %s %s' % (filename, dataname))

print_run_dict(run_dict)



this_time = datetime.now()
this_time_str = datetime.strftime(this_time, "%Y %m %d %H:%M:%S ")
print "Finish to build run_dir at %s" % (this_time_str)

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

#results_dict = {}
#for tuple in run_dict.keys():
#    results_dict[tuple] = analyse_properties(tuple)
#sys.exit()

pool = Pool()
results = pool.map(analyse_properties, run_dict.keys() )
results_dict = {}
filetuple = open('List-tuple.dat', 'w')
for tuple, result in zip( run_dict.keys(), results):
    filetuple.write('%s\n' % '  '.join(tuple))
    results_dict[tuple] = result
    properties = results_dict[tuple]
    for property in mean_properties:
        filename = create_file(property, title, properties[property + 'info'], 'Mean', tuple = tuple)
        os.system('mv %s %s' % (filename, dataname))
    for property in detail_properties:
        #print properties[property]
        properties[property] = average_dict(properties[property], properties['Number runs'])
        #print properties['Number runs']
        #print properties[property]

        filename = print_dict( properties[property], property, tuple  )
        os.system('mv %s %s' % (filename, dataname))
        #sys.exit()
filetuple.close()
os.system('mv List-tuple.dat %s' % dataname )







this_time = datetime.now()
this_time_str = datetime.strftime(this_time, "%Y %m %d %H:%M:%S ")
print "End of the analysis at %s" % (this_time_str)




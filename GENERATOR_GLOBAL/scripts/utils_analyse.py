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

def append_two_dict(dict1, dict2):
    result = {}
    if dict1 is None:
        result = dict2
    else:
        for key in dict1:
            value = dict1[key] + dict2[key]
            result[key] = value
    return result


def analyse_properties(tuple, run_dict, dict_properties):
    real_start_time = datetime.now()
    start_time = datetime.strftime( real_start_time, "%Y %m %d %H:%M:%S ")
    print "One worker for: %s starts at %s" % (tuple, start_time)
    list_dir = run_dict[tuple]
    properties_dict = {}
    for directory in list_dir:
        dir = FSSHRun(directory)
        for property in dict_properties.get('Detailled', []):
            prop = dir.extract(property)
            properties_dict[property] = sum_two_dict(properties_dict.get(property), prop)
        for property in dict_properties.get('Mean', []):
            prop = dir.extract(property)
            list = statistics( prop )
            if properties_dict.get(property) is None:
                properties_dict[property] = []
                properties_dict[property + 'info'] = ''
            properties_dict[property].append(list)
            line = '%s      %s\n' % ( directory, '    '.join(map(str, list )))
            properties_dict[property + 'info'] += line
        for property in dict_properties.get('Specific', []):
            prop = dir.extract(property)
            list = prop
            if properties_dict.get(property) is None:
                properties_dict[property] = []
                properties_dict[property + 'info'] = ''
            line = '%s      %s\n' % (directory, '    '.join(map(str, list)))
            properties_dict[property + 'info'] += line
        for property in dict_properties.get('Initial', []):
            prop = dir.extract(property, init = True)
            line = '%s      %s\n' % (directory, '    '.join(map(str, prop.values()[0])))
            if properties_dict.get(property) is None:
                properties_dict[property] = []
                properties_dict[property + 'info'] = ''
            properties_dict[property + 'info'] += line
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
    return filename


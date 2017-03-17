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


scripts = 'task234_parallel'

detail_properties = ['Surface populations', 'Adiabatic populations']
ratio_properties = ['Internal consistency ratio']
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
reorga= 0.300
free_energy = 0.00
nadiab = 2


dataname = 'data-%s-%s' % (scripts, title)
if not os.path.isdir(dataname):
    os.mkdir(dataname)
per_run_name = 'per-run-%s' % scripts
if not os.path.isdir(per_run_name):
    os.mkdir(per_run_name)


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

def analyse_properties(tuple):
    print "One core for: %s" % (tuple,)
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
    return properties_dict

#results_dict = {}
#for tuple in run_dict.keys():
#    results_dict[tuple] = analyse_properties(tuple)

pool = Pool()
results = pool.map(analyse_properties, run_dict.keys() )
results_dict = {}
for tuple, result in zip( run_dict.keys(), results):
    results_dict[tuple] = result

nfig =-1
for tuple in run_dict.keys():
    scaling = tuple[0]
    reversal = tuple[1]
    properties = results_dict[tuple]
    for property in mean_properties:
        filename = create_file(property, title, properties[property + 'info'], 'Mean')
        os.system('mv %s %s' % (filename, dataname))
    temperature = np.mean(np.array(properties['Temperature']).transpose()[0])
    coupling = np.mean(np.array(properties['Couplings']).transpose()[0]) * 27.211399  # coupling in eV

    boltzmann_ratio = calculate_boltzman_ratio(reorga, free_energy, coupling, temperature)

    nfig += 1
    plt.figure(nfig)
    for property in detail_properties:
        properties[property] = average_dict(properties[property], properties['Number runs'])
        plt.plot(sorted(properties[property]), [properties[property][time][0] for time in sorted(properties[property])],
                 label=property)
        if property == 'Surface populations':
            plt.plot(sorted(properties[property]),
                     [boltzmann_ratio for time in sorted(properties[property])], label='Boltzmann ratio')
    plt.xlabel('Time (fs)')
    plt.ylim([0, 1])
    plt.ylabel('Populations')
    plt.title('Populations vs time for reversal: %s and scaling value C = %s Ha' % (reversal, scaling))
    plt.legend()
    # plt.show()
    plt.savefig('%s/population_vs_time_%s_%s_%s.png' % (dataname, reversal, scaling, title))
    # sys.exit()
    # for time in sorted(properties[property]):
    # print time, properties[property][time]

    nfig += 1
    plt.figure(nfig)
    for property in ratio_properties:
        properties[property] = {}
        for time in sorted(properties['Surface populations']):
            if properties['Adiabatic populations'].get(time) is not None:
                properties[property][time] = \
                properties['Surface populations'][time] / \
                properties['Adiabatic populations'][time]
                # print time, properties[property][time], properties['Surface populations'][time], properties['Adiabatic populations'][time]
        plt.plot(sorted(properties[property]), [properties[property][time][0] for time in sorted(properties[property])],
                 label=property)
    plt.xlabel('Time (fs)')
    # plt.ylim([0, 1])
    plt.ylabel('Populations')
    plt.title('Ratio vs time for reversal: %s and scaling value C = %s Ha' % (reversal, scaling))
    plt.legend()
    plt.savefig('%s/ratio_vs_time_%s_%s_%s.png' % (dataname, reversal, scaling, title))
    # plt.show()
    # sys.exit()

print "End of the analysis"



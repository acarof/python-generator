import string, re, struct, sys, math, os, time
import numpy as np
import importlib, imp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from operator import itemgetter
from utils_scripts import *
from multiprocessing import Pool


scripts = 'task234_reversal_loop'

detail_properties = ['Surface populations', 'Adiabatic populations']
#histo_properties = ['Delta_E', 'Couplings', 'Populations']
histo_properties = ['Delta_E']
mean_properties = ['Temperature']
#specific_properties = ['FSSH']
specific_properties = []
total_properties = detail_properties + mean_properties + specific_properties + histo_properties
test = True

if test:
    dirlist = ['run-%d' % i for i in range(40)]
    title = 'TEST'
else:
    dirlist = os.listdir('.')
    name_bucket = os.getcwd().split('/')[-1]
    short_time = time.strftime("%y%m%d%H%M", time.localtime())
    title = '%s-%s' % (name_bucket, short_time,)


bin = 2
bin_histo = 50
#reorganization = 0.300
#free_energy = 0.00

properties_dict  = {}


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


nadiab = 2
os.system('cd ..')

def analyse_properties(directory):
#    if i % 100 == 0:
#        print time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    if 'run-' in directory and 'per' not in directory:
        print "Do %s" % directory
        os.chdir(directory)
        dir = FSSHRun(directory)
        [ scaling, reversal ] = dir.get_input_key(['SCALING_FACTOR','METHOD_REVERSAL'])
        if properties_dict.get(reversal) is None:
            properties_dict[reversal] = {}
        if properties_dict[reversal].get(scaling) is None:
            properties_dict[reversal][scaling] = {}
        property = 'Number runs'
        if properties_dict[reversal][scaling].get(property) is None:
            properties_dict[reversal][scaling][property] = 1
        else:
            properties_dict[reversal][scaling][property] += 1
        for property in (total_properties):
            prop = dir.extract(property)
            if property in detail_properties:
                properties_dict[reversal][scaling][property] = sum_two_dict(properties_dict[reversal][scaling].get(property), prop)
                #properties_.get(property).dict.get(density).append(prop)
                #dir.detailed_print(prop, property)

#            elif property in mean_properties:
#                list = statistics( prop )
#                properties_.get(property).dict.get(density).append(list)
#                line = '%s      %s\n' % ( directory, '    '.join(map(str, list )))
#                properties_.get(property).info += line
#            elif property in specific_properties:
#                list = prop
#                line = '%s      %s\n' % (directory, '    '.join(map(str, list)))
#                properties_.get(property).info += line
#            elif property in histo_properties:
#                #print property
#                #histo = histogram(prop)
#                properties_.get(property).dict.get(density).append(prop)
#                dir.detailed_print(prop, property)
#                #dir.histo_print(histo, property)
#        os.system('mv *run*.dat ../per-run')
        os.chdir('..')


pool = Pool()
pool.map(analyse_properties, dirlist)


def average_dict(dict1, number):
    result = {}
    for key in dict1:
        result[key] = np.array( dict1[key] ) / number
    return result

nfig =-1
ratio_properties = ['Internal consistency ratio']
for reversal in properties_dict:
    for scaling in properties_dict[reversal]:
        properties= properties_dict[reversal][scaling]
        nfig += 1
        plt.figure(nfig)
        for property in detail_properties:
            properties[property] = average_dict( properties[property], properties['Number runs'])
            plt.plot( sorted(properties[property]), [ properties[property][time][0] for time in  sorted(properties[property])], label = property)
        plt.xlabel('Time (fs)')
        plt.ylim([0,1])
        plt.ylabel('Populations')
        plt.title('Populations vs time for reversal: %s and scaling value C = %s Ha' % (reversal, scaling) )
        plt.legend()
        #plt.show()
        plt.savefig('%s/population_vs_time_%s_%s_%s.png' % (dataname, reversal, scaling, title))
        #sys.exit()
            #for time in sorted(properties_dict[reversal][scaling][property]):
                #print time, properties_dict[reversal][scaling][property][time]

        nfig += 1
        plt.figure(nfig)
        for property in ratio_properties:
            properties[property] = {}
            for time in sorted( properties_dict[reversal][scaling]['Surface populations'] ):
                if properties_dict[reversal][scaling]['Adiabatic populations'].get(time) is not None:
                    properties_dict[reversal][scaling][property][time] = properties_dict[reversal][scaling]['Surface populations'][time] / properties_dict[reversal][scaling]['Adiabatic populations'][time]
                 #print time, properties_dict[reversal][scaling][property][time], properties_dict[reversal][scaling]['Surface populations'][time], properties_dict[reversal][scaling]['Adiabatic populations'][time]
            plt.plot( sorted(properties[property]), [ properties[property][time][0] for time in  sorted(properties[property])], label = property)
        plt.xlabel('Time (fs)')
        #plt.ylim([0, 1])
        plt.ylabel('Populations')
        plt.title('Ratio vs time for reversal: %s and scaling value C = %s Ha' % (reversal, scaling))
        plt.legend()
        plt.savefig('%s/ratio_vs_time_%s_%s_%s.png' % (dataname, reversal, scaling, title))
        #plt.show()
        #sys.exit()


sys.exit()

for property in (mean_properties):
    create_file(property, title, properties_.get(property).info, 'Mean')
for property in (specific_properties):
    create_file(property, title, properties_.get(property).info, 'Spec')


print "End of the analysis"
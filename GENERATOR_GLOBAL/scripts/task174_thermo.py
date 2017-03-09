import string, re, struct, sys, math, os, time
import numpy as np
import importlib, imp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from operator import itemgetter
from utils_scripts import *



#densities = [0.01, 0.02, 0.03, 0.04]
#detail_properties = ['State']
detail_properties = []
#histo_properties = ['Delta_E', 'Couplings', 'Populations']
histo_properties = ['Delta_E']
mean_properties = ['Temperature']
#specific_properties = ['FSSH']
specific_properties = []
total_properties = detail_properties + mean_properties + specific_properties + histo_properties
test = False

if test:
    dirlist = ['run-%d' % i for i in range(40)]
else:
    dirlist = os.listdir('.')

bin = 2
bin_histo = 50
#reorganization = 0.300
#free_energy = 0.00

dict_natoms = {
    0.00: 12,
    0.0001: 19,
    0.0005: 75,
    0.001: 136,
    0.005: 1010,
    0.01: 218,
    0.02: 506,
    0.03: 712,
    0.04: 979,
    0.06: 317,
    0.08: 462,
    0.10: 651
}
dict_densities =  {v: k for k, v in dict_natoms.iteritems()}

densities = dict_densities.values()
my_densities = []


class PropDict(object):
    def __init__(self, property):
        self.prop = property
        self.dict = {}
        self.info = ''
        for density in densities:
            self.dict.update( { density : [] } )

properties_  = {}
for property in (total_properties) :
    properties_.update( {property :  PropDict(property)})

name_bucket = os.getcwd().split('/')[-1]
short_time = time.strftime("%y%m%d%H%M", time.localtime())
title = '%s-%s' % (name_bucket, short_time, )

dataname = 'data-%s' % title
if not os.path.isdir(dataname):
    os.mkdir(dataname)

if not os.path.isdir('per-run'):
    os.mkdir('per-run')


nadiab = 2
os.system('cd ..')
for i, directory in enumerate(dirlist):
    if i % 100 == 0:
        print time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    if 'run-' in directory:
        print "Do %s" % directory
        os.chdir(directory)
        dir = FSSHRun(directory)
        natoms = dir.extract('Atoms Number')
        density = dict_densities.get(natoms)
        if density not in my_densities:
            my_densities.append(density)
        for property in (total_properties):
            prop = dir.extract(property)
            if property in detail_properties:
                properties_.get(property).dict.get(density).append(prop)
                dir.detailed_print(prop, property)
            elif property in mean_properties:
                list = statistics( prop )
                properties_.get(property).dict.get(density).append(list)
                line = '%s      %s\n' % ( directory, '    '.join(map(str, list )))
                properties_.get(property).info += line
            elif property in specific_properties:
                list = prop
                line = '%s      %s\n' % (directory, '    '.join(map(str, list)))
                properties_.get(property).info += line
            elif property in histo_properties:
                #print property
                #histo = histogram(prop)
                properties_.get(property).dict.get(density).append(prop)
                dir.detailed_print(prop, property)
                #dir.histo_print(histo, property)
        os.system('mv *run*.dat ../per-run')
        os.chdir('..')

for property in (mean_properties):
    create_file(property, title, properties_.get(property).info, 'Mean')
for property in (specific_properties):
    create_file(property, title, properties_.get(property).info, 'Spec')



my_densities.sort()
for density in my_densities:
    print "Start density: %f" % density
    for property in histo_properties:
        bins, values, errors = histo_print(properties_.get(property).dict.get(density), property, title, bin, bin_histo, graph = True)
        fig = 1
        plt.figure(fig, figsize=(6, 5))
        plt.plot(bins, values, drawstyle='steps', label=r'%s (atoms/$\AA^3$)' % density)
        plt.title('Distribution of site energies')
        plt.xlabel(r'$\Delta E$ (Ha)')
        plt.ylabel('Density of probability (Ha$^{-1})$')
        plt.yscale('log')
        #plt.ylabel(r'Rate (fs$^{-1}$)')
        #plt.legend()
        #plt.legend(bbox_to_anchor=(1.6, 1))

        temperatures = np.array(properties_.get('Temperature').dict.get(density)).transpose()[0]
        mean_temp = np.mean(temperatures)
        free_energy = -  mean_temp * np.log(values)

        fig = 2
        plt.figure(fig, figsize=(6, 5))
        plt.plot( bins, free_energy, drawstyle='steps', label=r'%s (atoms/$\AA^3$)' % density)
        plt.title('Distribution of site energies')
        plt.xlabel(r'$\Delta E$ (Ha)')
        plt.ylabel('Landau energy (K)')


#plt.savefig('histogram_%s_%s.png' % (property, title), bbox_inches='tight')



os.system('mv *.png %s' % dataname)

plt.show()
#plt.close()



print "End of the analysis"
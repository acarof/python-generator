import string, re, struct, sys, math, os, time
import numpy as np
#import importlib, imp
#import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit

def read_xyz_file(filename):
    file = open(filename, 'r')
    results = {}
    for line in file.readlines():
        if 'nadiab' in line:
            pass
        elif 'time' in line:
            time =  float(line.split()[5])
            results.update({time : []})
        else:
            results.get(time).append([ float(x) for x in line.split()])
    file.close()
    return results


def extract_coupling(filename = 'run-hamilt-1.xyz'):
    hamiltonians = read_xyz_file(filename)

    couplings = []
    for time in hamiltonians:
        couplings.append(hamiltonians.get(time)[0][3])
    return couplings


os.system('cd ..')
print "Coupling in meV"
means = []
fileout = open('electronic_coupling.dat', 'w')
for directory in os.listdir('.'):
    if 'run' in directory:
        os.chdir(directory)
        couplings = extract_coupling()
        mean = np.sqrt( np.mean( np.square(couplings) ) )*27.211399*1000
        means.append(np.sqrt( mean ) )
        print "%s %s" % (directory, mean)
        fileout.write("%s %s\n" % (directory, mean))

        os.chdir('..')

fileout.close()
print "Average coupling is: %s meV" % np.mean(means)
print "End of the analysis"
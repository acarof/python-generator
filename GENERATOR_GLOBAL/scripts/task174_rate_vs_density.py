import string, re, struct, sys, math, os, time
import numpy as np
import importlib, imp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit




def extract_population(filename = 'run-coeff-1.xyz'):
    file = open(filename,'r')
    line = file.readline()
    nadiab = int(line.split()[3].strip(';'))
    norbital = int(line.split()[6])

    pop = [[] for x in range(nadiab*norbital)]
    times = []
    for line in file.readlines():
        if 'nadiab' in line:
            pass
        elif 'time' in line:
            times.append(float(line.split()[5]))
        else:
            index = norbital * (int(line.split()[0]) - 1) + int(line.split()[1]) - 1
            coefficients = np.complex(float(line.split()[2]), float(line.split()[3]))
            population = np.square(np.absolute(coefficients)       )
            pop[index].append(population)
    return times, pop, nadiab*norbital


def extract_phase(filename = 'run-coeff-1.xyz'):
    file = open(filename,'r')
    line = file.readline()
    nadiab = int(line.split()[3].strip(';'))
    norbital = int(line.split()[6])

    quant = [[] for x in range(nadiab*norbital)]
    times = []
    for line in file.readlines():
        if 'nadiab' in line:
            pass
        elif 'time' in line:
            times.append(float(line.split()[5]))
        else:
            index = norbital * (int(line.split()[0]) - 1) + int(line.split()[1]) - 1
            coefficients = np.complex(float(line.split()[2]), float(line.split()[3]))
            quantity = np.angle(coefficients)
            quant[index].append(quantity)
    return times, quant, nadiab*norbital

def func(x, a, b, c):
    return a*np.exp(- b*x) + c


def get_atoms_number(filename = 'run.inp'):
    file = open(filename, 'r')
    for line in file.readlines():
        if  re.search('NUMBER_OF_ATOMS', line):
            natoms = int(line.split()[1])
    return natoms

list_density = [0.01, 0.02, 0.03, 0.04]
density_dict = {
    218 : 0.01,
    506 : 0.02,
    712 : 0.03,
    979 : 0.04
}
stimes_dict = {
    0.01: [],
    0.02: [],
    0.03: [],
    0.04: []
}
spop_dict = {
    0.01: [],
    0.02: [],
    0.03: [],
    0.04: []
}


os.system('cd ..')
for directory in os.listdir('.'):
    if 'run' in directory:
        os.chdir(directory)
        natoms = get_atoms_number()
        density = density_dict.get(natoms)
        print density
        times, pop, nadiab = extract_population()
        angle = extract_phase()[1]
        stimes_dict.get(density).append(times)
        spop_dict.get(density).append(pop)
        os.chdir('..')


for density in list_density:
    stimes = np.array(stimes_dict.get(density))
    spop   = np.array(spop_dict.get(density))
    meantimes = np.mean(stimes, axis=0)
    meanpop = np.mean(spop, axis = 0)
    plt.figure(1)
    for adiab in range(nadiab-1):
        plt.plot(meantimes, meanpop[adiab], label = adiab)
        #popt, pcov = curve_fit(func, meantimes, meanpop[adiab])
        #plt.plot(meantimes, func(meantimes, *popt), label = 'Fitted curve')
        #print "The decay rate for state %d is %f fs-1" %(adiab, popt[1])
    plt.plot(meantimes, np.sum(meanpop, axis = 0))
plt.title('Population evolution')
plt.xlabel('Time (fs)')
plt.locator_params(axis='x', nbins=5)
plt.ylabel('Population')
#plt.ylim([5, 36])
plt.legend(bbox_to_anchor=(2.01, 1))
plt.show()
plt.clf()
plt.cla()
plt.close()


print "End of the analysis"
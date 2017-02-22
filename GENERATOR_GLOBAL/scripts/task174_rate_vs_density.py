import string, re, struct, sys, math, os, time
import numpy as np
import importlib, imp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from operator import itemgetter

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


def mean_over_run(list_of_dict):
    means = []
    for time in list_of_dict[0]:
        mean = [0] * len( list_of_dict[0].get(time) )
        counter = 0
        for dict_ in list_of_dict:
            mean = mean + np.array( dict_.get(time) )
            counter = counter + 1
        mean = mean / counter
        mean = np.insert(mean, 0, time)
        means.append(mean)
    means = sorted(means, key=itemgetter(0))
    means = map(list, zip(*means))
    return means




def extract_coupling(filename = 'run-hamilt-1.xyz'):
    hamiltonians = read_xyz_file(filename)

    couplings = []
    for time in hamiltonians:
        couplings.append(hamiltonians.get(time)[0][3])
    return couplings

def extract_population(filename = 'run-coeff-1.xyz'):
    coefficients = read_xyz_file(filename)

    populations = {}
    for time in coefficients:
        list = coefficients.get(time)
        pop = []
        for i, coeff in enumerate(list):
            pop.append(np.square( np.absolute( np.complex( coeff[2], coeff[3]) ) ) )
        populations.update({time : pop})
    return populations


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
scoupling_dict = {
    0.01: [],
    0.02: [],
    0.03: [],
    0.04: []
}


nadiab = 2
os.system('cd ..')
for directory in os.listdir('.'):
    if 'run' in directory:
        print "Do %s" % directory
        os.chdir(directory)
        natoms = get_atoms_number()
        density = density_dict.get(natoms)
        pop = extract_population()
        coupling = extract_coupling()
        mean_coupling = np.sqrt(np.mean(np.square(coupling))) * 27.211399 * 1000
        scoupling_dict.get(density).append(mean_coupling)
        spop_dict.get(density).append(pop)
        os.chdir('..')




def calculate_barrier(reorganization, free_energy):
    energy_barrier = np.square( ( reorganization + free_energy ) ) / ( 4 * reorganization)
    return  energy_barrier


def calculate_marcus_na_rate(coupling, reorganization, free_energy, temperature):
    #Convert everything in atomic units
    coupling = coupling / 27.211399
    reorganization = reorganization / 27.211399
    free_energy = free_energy /     27.211399
    temperature = temperature / 315777.09

    energy_barrier = calculate_barrier(reorganization, free_energy)

    rate = (2 * np.pi) * np.square(coupling) * np.exp( - energy_barrier / temperature) / np.sqrt( 4 * np.pi * reorganization * temperature)
    rate = rate / 0.024188
    return rate



def func(x, a, b, c):
    return a*np.exp(- b*x) + c

reorganization = 0.300
free_energy = 0.00
temperature = 300

plt.figure(1)
rates=[]
marcus_rates = []
for density in list_density:
    spop   = np.array(spop_dict.get(density))
    scoupling = np.array(scoupling_dict.get(density))
    pop = mean_over_run(spop)
    for adiab in range(1, nadiab):
        plt.plot(pop[0], pop[adiab], label = adiab)
        try:
            popt, pcov = curve_fit(func, pop[0], pop[adiab])
        except:
            popt = 0.0, 0.0, 0.0
        rate = popt[1] / 2
        rates.append( rate)
        plt.plot( pop[0], [ func( x, *popt) for x in pop[0] ], label = 'Fitted curve')
        print "The decay rate for state %d is %f fs-1" %(adiab, rate)
    mean_coupling =  np.mean( scoupling )
    print "The average coupling is %f meV " %  mean_coupling
    marcus_rate = calculate_marcus_na_rate(mean_coupling / 1000, reorganization, free_energy, temperature)
    marcus_rates.append( marcus_rate)
    print "The Marcus rate is  %f fs-1" % marcus_rate
plt.title('Population evolution')
plt.xlabel('Time (fs)')
plt.locator_params(axis='x', nbins=5)
plt.ylabel('Population')
#plt.ylim([5, 36])
plt.legend(bbox_to_anchor=(2.01, 1))






plt.figure(2)
plt.plot(list_density, rates)
plt.plot(list_density, marcus_rates)
plt.xlabel(r'Density (atoms/$\AA^3$)')
plt.yscale('log')
plt.ylabel(r'Rate (fs$^{-1}$)')

plt.show()
plt.close()

print "End of the analysis"
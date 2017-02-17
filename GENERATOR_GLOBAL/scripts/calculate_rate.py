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


os.system('cd ..')
stimes = []
spop = []
sangle=[]
for directory in os.listdir('.'):
    if 'run' in directory:
        os.chdir(directory)
        times, pop, nadiab = extract_population()
        angle = extract_phase()[1]
        sangle.append(angle)
        stimes.append(times)
        spop.append(pop)
        os.chdir('..')

meantimes = np.mean(stimes, axis=0)

meanpop = np.mean(spop, axis = 0)
plt.figure(1)
for adiab in range(nadiab-1):
    plt.plot(meantimes, meanpop[adiab], label = adiab)
    popt, pcov = curve_fit(func, meantimes, meanpop[adiab])
    plt.plot(meantimes, func(meantimes, *popt), 'r-', label = 'Fitted curve')
    print "The decay rate for state %d is %f fs-1" %(adiab, popt[1])
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


plt.figure(2)
meanangle = np.mean(sangle, axis = 0)
for adiab in range(1, nadiab):
    plt.plot(meantimes, meanangle[adiab] - meanangle[0], label = adiab)
plt.title('Phase evolution')
plt.xlabel('Time (fs)')
plt.locator_params(axis='x', nbins=5)
plt.ylabel('Phase (rad)')
#plt.ylim([5, 36])
plt.legend(bbox_to_anchor=(2.01, 1))
plt.show()
plt.clf()
plt.cla()
plt.close()

print "End of the analysis"
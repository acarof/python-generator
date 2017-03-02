import string, re, struct, sys, math, os, time
import numpy as np
from operator import itemgetter

units = {
    'Couplings' : 'Ha',
    'Populations' : '',
    'Temperature' : 'K',
    'Delta_E'     : 'Ha'
}

headers ={
    'FSSH' : '# Attempts  Kin-passed Successful Decoherence\n '
}


class FSSHRun(object):
    def __init__(self, name):
        self.name = name

    def extract(self, property):
        if ( property == 'Atoms Number'):
            return self._get_atoms_number()
        if ( property == 'Couplings'):
            return self._extract_coupling()
        if ( property == 'Populations'):
            return self._extract_population()
        if ( property == 'Temperature'):
            return self._extract_energies('Temperature')
        if ( property == 'FSSH'):
            return self._extract_fssh()
        if ( property == 'Delta_E'):
            return self._extract_delta_e()


    def detailed_print(self, prop, propname):
        filename = '%s-%s.dat' % (propname, self.name)
        file = open(filename, 'w')
        file.write( '# Time  %s(%s)\n' % (propname, units.get(propname)) )
        for time in sorted(prop):
            line = '%20.10f ' % time
            for element in prop.get(time):
                line = line + ('%20.10f ' % element)
            file.write(line + '\n')
        file.close()

    def histo_print(self, histo, propname):
        filename = '%s-histo-%s.dat' % (propname, self.name)
        file = open(filename, 'w')
        file.write('# %s(%s) Density of probability\n' % (propname, units.get(propname)) )
        for value, bin in zip(histo[0], histo[1]):
            line = '%20.10f %20.10f ' % (bin, value)
            file.write(line + '\n')
        file.close()


    def _extract_fssh(self, filename = 'run-sh-1.log'):
        results = []
        for property in ['ATTEMPT', 'PASSED', 'SUCCESS','DECOHERENCE']:
            with open(filename) as f:
                contents = f.read()
                count = contents.count(property)
            results.append(count)
        return results


    def _get_atoms_number(self, filename = 'run.inp'):
        file = open(filename, 'r')
        for line in file.readlines():
            if  re.search('NUMBER_OF_ATOMS', line):
                natoms = int(line.split()[1])
        return natoms


    def _read_xyz_file(self, filename):
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

    def _read_ener_file(self):
        file = open('run-1.ener', 'r')
        results = {}
        for line in file.readlines():
            if '#' in line:
                pass
            else:
                time = float(line.split()[1])
                list = map(float, line.split() )
                results.update({time : [] })
                results[time] = results.get(time) + list
        file.close()
        return results

    def _extract_energies(self, property):
        properties = {
            'Steps' : 0,
            'Time'  : 1,
            'Kinetic' : 2,
            'Temperature' : 3,
            'Potential'   : 4,
            'Conserved'   : 5,
            'UsedTimes'   : 6

        }
        column = properties.get(property)
        if column is None:
            raise NotImplementedError
        energies = self._read_ener_file()
        prop = {}
        for time in energies:
            prop[time] = energies.get(time)[column]
            #prop.append( energies.get(time)[column] )
        return prop

    def _extract_delta_e(self, filename = 'run-hamilt-1.xyz'):
        hamiltonians = self._read_xyz_file(filename)
        delta_e = {}
        for time in hamiltonians:
            delta_e[time] = []
            for i in range(1, len(hamiltonians.get(time))):
                delta_e[time].append(
                    hamiltonians.get(time)[i][i+2] - hamiltonians.get(time)[i-1][i-1+2]
                )
        return delta_e

    def _extract_coupling(self, filename = 'run-hamilt-1.xyz'):
        hamiltonians = self._read_xyz_file(filename)
        couplings = {}
        for time in hamiltonians:
            couplings[time] = hamiltonians.get(time)[0][3]
            #couplings.append(hamiltonians.get(time)[0][3])
        return couplings

    def _extract_population(self, filename = 'run-coeff-1.xyz'):
        coefficients = self._read_xyz_file(filename)
        populations = {}
        for time in coefficients:
            list = coefficients.get(time)
            pop = []
            for i, coeff in enumerate(list):
                pop.append(np.square( np.absolute( np.complex( coeff[2], coeff[3]) ) ) )
            populations.update({time : pop})
        return populations




def statistics(prop):
    times = []
    props = []
    for time in sorted(prop):
        times.append(time)
        props.append(prop.get(time))
    mean = np.mean(props)
    std = np.std(props)
    drift = np.polyfit(times, props, 1)[0]
    qmean = np.sqrt(np.mean(np.square(props)))
    return [mean, std, drift, qmean]

def histogram(prop):
    props = []
    for time in prop:
        props.append(prop.get(time))
    return np.histogram(props, normed=True)


def mean_over_run(list_of_dict):
    means = []
    for time in list_of_dict[0]:
        mean = [0] * len(list_of_dict[0].get(time))
        counter = 0
        for dict_ in list_of_dict:
            mean = mean + np.array(dict_.get(time))
            counter = counter + 1
        mean = mean / counter
        mean = np.insert(mean, 0, time)
        means.append(mean)
    means = sorted(means, key=itemgetter(0))
    means = map(list, zip(*means))
    return means


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

def expo(x, a, b, c):
    return a*np.exp(- b*x) + c

def create_file(property, bucket, text, label = 'None'):
    filename = '%s-%s.dat' % (property, bucket)
    file = open(filename, 'w')
    unit = units.get(property)
    if (label == 'Mean'):
        replace = (unit,) * 4
        header = '# Run  Mean (%s)   Drift (%s/fs)  Std (%s)  QMean (%s) \n' %  replace
    elif (label == 'Spec'):
        header = headers.get(property)
    else:
        header = None
    if header is not None:
        file.write(header)
    file.write(text)
    file.close()
    os.system('mv %s data/' % filename)
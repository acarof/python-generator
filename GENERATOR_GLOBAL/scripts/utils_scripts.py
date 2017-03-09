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
        if ( property == 'State'):
            return self._extract_state()


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


    def _extract_state(self, filename = 'run-sh-1.log' ):
        file = open(filename, 'r')
        results = {}
        for line in file.readlines():
            if 'time' in line:
                time = float(line.split()[5])
                results.update({time : []})
            elif 'Final' in line:
                state = int( line.split()[3] )
                results.get(time).append(state)
        file.close()
        return results



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
            prop[time] = [energies.get(time)[column]]
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
            couplings[time] = [hamiltonians.get(time)[0][3]]
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
    drift = np.polyfit(times, props, 1)[0][0]
    qmean = np.sqrt(np.mean(np.square(props)))
    return [mean, std, drift, qmean]

def histogram(prop):
    props = []
    for time in prop:
        props.append(prop.get(time))
    return np.histogram(props, normed=True)

def histo_print(list_of_dict, propname, title, bin = 1, bin_histo = 10, list_state = None, nadiab = 2, graph = False):
    if list_state is not None:
        histo_print_per_state(list_of_dict, propname, title, bin, bin_histo, list_state, nadiab)
    else:
        filename = '%s-histo-%s.dat' % (propname, title)
        file = open(filename, 'w')
        file.write('# %s(%s) Density_of_probability  Error\n' % (propname, units.get(propname)))
        n = int( len(list_of_dict) / bin )
        llod = [list_of_dict[i:i + n ] for i in xrange(0, len(list_of_dict), n)]
        slist = [[] for x in xrange(bin)]
        for i in range(bin):
            slist[i] = []
            for dict in llod[i]:
                for value in dict.values():
                    slist[i] += value
        max_ =  max([sublist[-1] for sublist in slist])
        min_ =  min([sublist[-1] for sublist in slist])
        histos = []
        bins = []
        for i in range(bin):
            bins.append( np.histogram(slist[i], bins = bin_histo, range=(min_, max_), normed = True) [1] )
            histos.append( np.histogram(slist[i], bins = bin_histo, range=(min_, max_), normed = True) [0] )
        binhs = []
        values = []
        errors = []
        for binh, value, error in zip( np.mean( bins, axis = 0),
                                       np.append( np.mean( histos, axis = 0), [0]),
                                       np.append( np.std( histos, axis=0, ddof=1 ), [0] )):
            line = '%20.10f %20.10f %20.10f ' % (binh, value, error)
            binhs.append(binh)
            values.append(value)
            errors.append(error)
            file.write(line + '\n')
        file.close()
        os.system('mv %s data-%s' %  (filename, title))
        if graph:
            return  binhs, values, errors


def histo_print_per_state(list_of_dict, propname, title, bin = 1, bin_histo = 10, list_state = None, nadiab = 2):
    n = int( len(list_of_dict) / bin )
    llod = [list_of_dict[i:i + n ] for i in xrange(0, len(list_of_dict), n)]
    llos = [list_state[i:i + n ] for i in xrange(0, len(list_state), n)]#

    states_dict = { i :  [[] for x in xrange(bin)] for i in range(1, nadiab+1)}

    for i in range(bin):
        for dict, states in zip(llod[i], llos[i]):
            for time in dict:
                if states.get(time) is not None:
                    states_dict[  states.get(time)[0] ][i] += dict.get(time)

    for state in range(1, nadiab+1):
        filename = '%s-histo-state-%d-%s.dat' % (propname, state, title)
        file = open(filename, 'w')
        file.write('# %s(%s) Density_of_probability  Error\n' % (propname, units.get(propname)))
        slist = states_dict[ state ]
        maxlist = []
        minlist = []
        for sublist in slist:
            try:
                maxlist.append(max(sublist))
                minlist.append(min(sublist))
            except:
                pass
        max_ = max(maxlist)
        min_ = min(minlist)
        histos = []
        bins = []
        for i in range(bin):
            if slist[i]:
                bins.append( np.histogram(slist[i], bins = bin_histo, range=(min_, max_), normed = True) [1] )
                histos.append( np.histogram(slist[i], bins = bin_histo, range=(min_, max_), normed = True) [0] )
        for binh, value, error in zip( np.mean( bins, axis = 0),
                                       np.append( np.mean( histos, axis = 0), [0]),
                                       np.append( np.std( histos, axis=0, ddof=1 ), [0] )):
            line = '%20.10f %20.10f %20.10f ' % (binh, value, error)
            file.write(line + '\n')
        file.close()
        os.system('mv %s data-%s' %  (filename, title))




def mean_over_run(list_of_dict, bin =1):
    n = int( len(list_of_dict) / bin )
    llod = [list_of_dict[i:i + n ] for i in xrange(0, len(list_of_dict), n)]
    means = [[] for x in xrange(bin)]
    for i in range(bin):
        submeans = []
        for time in llod[i][0]:
            mean = [0] * len(llod[i][0].get(time))
            counter = 0
            for dict_ in llod[i]:
                mean = mean + np.array(dict_.get(time))
                counter = counter + 1
            mean = mean / counter
            mean = np.insert(mean, 0, time)
            submeans.append(mean)
        submeans = sorted(submeans, key=itemgetter(0))
        submeans = map(list, zip(*submeans))
        means[i] = submeans
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

def expo_constraint(x, b):
    return 0.5*np.exp(- b*x) + 0.5

def expo_free(x, b, a, c):
    return a*np.exp(- b*x) + c

def log_constraint(x, b):
    return -b*x

def log_free(x, b, a, c):
    return np.log( a*np.exp(- b*x) + c )

def create_file(property, title, text, label = 'None'):
    filename = '%s-%s.dat' % (property, title)
    file = open(filename, 'w')
    unit = units.get(property)
    if (label == 'Mean'):
        replace = (unit,) * 4
        header = '# Run  Mean (%s)  Std (%s)  Drift (%s/fs)   QMean (%s) \n' %  replace
    elif (label == 'Spec'):
        header = headers.get(property)
    else:
        header = None
    if header is not None:
        file.write(header)
    file.write(text)
    file.close()
    os.system('mv %s data-%s' % (filename, title))

def print_dat(array, label, title):
    new = map(list, zip(*array))
    filename = label + '.dat'
    file = open(filename, 'w')
    for line in new:
        file.write('    '.join(map(str, line )) + '\n')
    file.close()
    os.system('mv %s data-%s' % (filename, title))



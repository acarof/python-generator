import string, re, struct, sys, math, os, time
import numpy as np
from operator import itemgetter

units = {
    'Couplings': 'Ha',
    'Populations': '',
    'Temperature': 'K',
    'Delta_E': 'Ha'
}

headers = {
    'FSSH': '# Attempts  Kin-passed Successful Decoherence\n ',
    'Detailed-FSSH' : '# Attempts  Kin-passed Successful Decoherence Ups Downs\n'
}


class FSSHRun(object):
    def __init__(self, name):
        self.name = name

    def extract(self, property, init=False, **kwargs):
        if init:
            self.init = True
            self.timestep  = float(self.get_input_key(['\TIMESTEP'])[0])
        else:
            self.init = False
        if (property == 'Atoms Number'):
            return self._get_atoms_number()
        elif (property == 'Couplings'):
            return self._extract_coupling()
        elif (property == 'Populations'):
            return self._extract_population()
        elif (property == 'Temperature'):
            return self._extract_energies('Temperature')
        elif (property == 'Total-energy'):
            return self._extract_energies('Conserved')
        elif (property == 'FSSH'):
            return self._extract_fssh()
        elif property == 'Detailed-FSSH':
            return self._extract_detailed_fssh()
        elif (property == 'Delta_E'):
            return self._extract_delta_e()
        elif (property == 'State'):
            return self._extract_state()
        elif (property == 'Forces'):
            return self._extract_forces(**kwargs)
        elif (property == 'NACE'):
            return self._extract_coupling(filename='run-nace-1.xyz')
        elif (property == 'NACV'):
            return self._extract_nacv()
        elif property == 'Surface-populations':
            return self._extract_surface_population()
        elif 'Surface-populations-' in property:
            surface = int(property.replace("Surface-populations-", ""))
            return  self._extract_surface_population(surface=surface)
        elif property == 'Adiabatic-populations':
            return self._extract_adiabatic_population()
        elif property == 'Adiabatic-energies':
            return self._extract_adiabatic_energies()
        elif property == 'MSD':
            return self._extract_msd()
        else:
            print "Extraction of %s not implemented" % property
            raise SystemExit

    def _extract_msd(self):
        populations = self._extract_population()
        com = self._extract_com()
        result = {}
        for time, pop in populations.items():
            if com.get(time):
                msd = 0.0
                for (x,y) in zip(pop, com[time]):
                    msd += x*y*y
                result[time] = msd
        return result

    def _extract_com(self, filename = 'run-pos-1.xyz'):
        populations = self._extract_population()
        results = {}
        for time, pop in populations.items():
            results[time] = [ x*6.038 for x in range(len(pop))]
        return results


    def _extract_state(self, filename='run-sh-1.log'):
        file = open('%s/%s' % (self.name, filename), 'r')
        results = {}
        for line in file.readlines():
            if ', time' in line:
                time = float(line.split()[5])
                if self.init:
                    if time > self.timestep:
                        file.close()
                        return results
                results.update({time: []})
            elif 'Final' in line:
                state = int(line.split()[3])
                results.get(time).append(state)
        file.close()
        return results

    def _extract_surface_population(self, surface = None):
        def _state_to_surface_pop(state, nadiab):
            list = [0.00] * nadiab
            list[state - 1] = 1.00
            return list

        state = self._extract_state()
        results = {}
        nadiab = int( self.get_input_key(['NUMBER_DIABATIC_STATES'])[0] )
        for time in state:
            list_surface = _state_to_surface_pop(state[time][0], nadiab)
            if surface is None:
                results[time] = list_surface
            else:
                results[time] = [list_surface[surface]]
        return results

    def _extract_fssh(self, filename='run-sh-1.log'):
        results = []
        for property in ['ATTEMPT', 'PASSED', 'SUCCESS', 'DECOHERENCE!']:
            with open('%s/%s' % (self.name, filename)) as f:
                contents = f.read()
                count = contents.count(property)
            results.append(count)
        return results

    def _extract_detailed_fssh(self, filename='run-sh-1.log'):
        def _extract_deco_per_state(filename='run-sh-1.log'):
            file = open('%s/%s' % (self.name, filename), 'r')
            nadiab = int(self.get_input_key(['NUMBER_DIABATIC_STATES'])[0])
            results = [0] * nadiab
            for line in file.readlines():
                if 'Final' in line:
                    state = int(line.split()[3])
                elif 'DECOHERENCE!' in line:
                    results[state - 1] += 1
            file.close()
            return results

        def _extract_up_down_hop(filename='run-sh-1.log'):
            results = [0, 0]
            states = self._extract_state(filename)
            previous = states[sorted(states)[0]][0]
            for time in sorted(states)[1:]:
                present = states[time][0]
                if (present - previous) > 0:
                    results[0] += 1
                elif (present - previous) < 0:
                    results[1] += 1
                previous = present
            return results


        results = self._extract_fssh(filename)
        results += _extract_up_down_hop(filename)
        deco = self.get_input_key(['DECOHERENCE_CORRECTIONS'])
        if deco[0] == 'INSTANT_COLLAPSE':
            results += _extract_deco_per_state(filename)
        return results


    def get_input_key(self, key_list, filename='run.inp'):
        file = open('%s/%s' % (self.name, filename), 'r')
        values = []
        preval = {}
        for line in file.readlines():
            for key in key_list:
                if re.search(key, line):
                    preval[key] = (line.split()[-1])
        for key in key_list:
            values.append( preval[key] )
        return tuple(values)

    def _extract_nacv(self, filename='run-nacv-1.xyz'):
        file = open('%s/%s' % (self.name, filename), 'r')
        results = {}
        for line in file.readlines():
            if 'nadiab' in line or len(line.split()) == 1:
                pass
            elif 'time' in line:
                pretime = line.split()[5]
                if pretime[-1] == ',':
                    pretime = pretime[:-1]
                time = float(pretime)
                results.update({time: {}})
            elif 'Atom' in line:
                atom = int( line.split()[2] )
                if results.get(time).get(atom) is None:
                    results.get(time).update({atom: {}})
                dim = int( line.split()[5] )
                results.get(time).get(atom).update({dim: []})
            else:
                results.get(time).get(atom).get(dim).append( [ float(x) for x in line.split()] )
        file.close()
        return results

    def _read_xyz_file(self, filename):
        def _convert_float(s):
            try:
                return float(s)
            except ValueError:
                return None

        def _convert_int(s):
            try:
                return int(s)
            except ValueError:
                return None

        file = open('%s/%s' % (self.name, filename), 'r')
        results = {}
        for line in file.readlines():
            if 'nadiab' in line or len(line.split()) == 1:
                pass
            elif 'time' in line:
                pretime = line.split()[5]
                if pretime[-1] == ',':
                    pretime = pretime[:-1]
                time = float(pretime)
                if self.init:
                    if (time > self.timestep):
                        file.close()
                        return results
                results.update({time: []})
            else:
                list = []
                for x in line.split():
                    if _convert_float(x) is None:
                        list.append(x)
                    else:
                        if _convert_int(x) is None:
                            list.append(_convert_float(x))
                        else:
                            list.append(_convert_int(x))

                results.get(time).append(list)
        file.close()
        return results

    def _extract_adiabatic_population(self):
        hamiltonian = self._read_xyz_file('run-hamilt-1.xyz')
        diabat_coeff = self._read_xyz_file('run-coeff-1.xyz')
        result = {}
        for time in hamiltonian:
            if diabat_coeff.get(time) is not None:
                hamilt = np.transpose( hamiltonian[time])[2:]
                eigenvalues, eigenvectors = np.linalg.eig( hamilt )
                idx = eigenvalues.argsort()
                eigenvalues = eigenvalues[idx]
                eigenvectors = eigenvectors[:, idx]
                diabat = [ np.complex(coeff[2], coeff[3]) for coeff in diabat_coeff[time] ]
                adiabat = eigenvectors.transpose().dot(diabat)
                populations = [ np.absolute(x)**2 for x in adiabat]
                result[time] = populations
        return result

    def _extract_adiabatic_energies(self):
        hamiltonian = self._read_xyz_file('run-hamilt-1.xyz')
        result = {}
        for time in hamiltonian:
            hamilt = np.transpose( hamiltonian[time])[2:]
            eigenvalues, eigenvectors = np.linalg.eig( hamilt )
            idx = eigenvalues.argsort()
            eigenvalues = eigenvalues[idx]
            result[time] = eigenvalues
        return result

    def _extract_energies(self, property):
        def _read_ener_file():
            filename = 'run-1.ener'
            file = open('%s/%s' % (self.name, filename), 'r')
            results = {}
            for line in file.readlines():
                if '#' in line:
                    pass
                else:
                    time = float(line.split()[1])
                    list = map(float, line.split())
                    results.update({time: []})
                    results[time] = results.get(time) + list
            file.close()
            return results

        properties = {
            'Steps': 0,
            'Time': 1,
            'Kinetic': 2,
            'Temperature': 3,
            'Potential': 4,
            'Conserved': 5,
            'UsedTimes': 6

        }
        column = properties.get(property)
        if column is None:
            raise NotImplementedError
        energies = _read_ener_file()
        prop = {}
        for time in energies:
            prop[time] = [energies.get(time)[column]]
            # prop.append( energies.get(time)[column] )
        return prop

    def _extract_forces(self, filename='run-frc-1.xyz'):
        forces = self._read_xyz_file(filename)
        frc = {}
        for time in forces:
            frc[time] = []
            for i in range(len(forces.get(time))):
                frc[time].append(forces.get(time)[i][1:4])
        return frc

    def _extract_delta_e(self, filename='run-hamilt-1.xyz'):
        hamiltonians = self._read_xyz_file(filename)
        delta_e = {}
        for time in hamiltonians:
            delta_e[time] = []
            for i in range(1, len(hamiltonians.get(time))):
                delta_e[time].append(
                    hamiltonians.get(time)[i][i + 2] - hamiltonians.get(time)[i - 1][i - 1 + 2]
                )
        return delta_e

    def _extract_coupling(self, filename='run-hamilt-1.xyz'):
        hamiltonians = self._read_xyz_file(filename)
        couplings = {}
        for time in hamiltonians:
            couplings[time] = [hamiltonians.get(time)[0][3]]
            # couplings.append(hamiltonians.get(time)[0][3])
        return couplings

    def _extract_population(self, filename='run-coeff-1.xyz'):
        coefficients = self._read_xyz_file(filename)
        populations = {}
        for time in coefficients:
            list = coefficients.get(time)
            pop = []
            for i, coeff in enumerate(list):
                pop.append(np.square(np.absolute(np.complex(coeff[2], coeff[3]))))
            populations.update({time: pop})
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


def histo_print(list_of_dict, propname, title, bin=1, bin_histo=10, list_state=None, nadiab=2, graph=False):
    if list_state is not None:
        histo_print_per_state(list_of_dict, propname, title, bin, bin_histo, list_state, nadiab)
    else:
        filename = '%s-histo-%s.dat' % (propname, title)
        file = open(filename, 'w')
        file.write('# %s(%s) Density_of_probability  Error\n' % (propname, units.get(propname)))
        n = int(len(list_of_dict) / bin)
        llod = [list_of_dict[i:i + n] for i in xrange(0, len(list_of_dict), n)]
        slist = [[] for x in xrange(bin)]
        for i in range(bin):
            slist[i] = []
            for dict in llod[i]:
                for value in dict.values():
                    slist[i] += value
        max_ = max([sublist[-1] for sublist in slist])
        min_ = min([sublist[-1] for sublist in slist])
        histos = []
        bins = []
        for i in range(bin):
            bins.append(np.histogram(slist[i], bins=bin_histo, range=(min_, max_), normed=True)[1])
            histos.append(np.histogram(slist[i], bins=bin_histo, range=(min_, max_), normed=True)[0])
        binhs = []
        values = []
        errors = []
        for binh, value, error in zip(np.mean(bins, axis=0),
                                      np.append(np.mean(histos, axis=0), [0]),
                                      np.append(np.std(histos, axis=0, ddof=1), [0])):
            line = '%20.10f %20.10f %20.10f ' % (binh, value, error)
            binhs.append(binh)
            values.append(value)
            errors.append(error)
            file.write(line + '\n')
        file.close()
        os.system('mv %s data-%s' % (filename, title))
        if graph:
            return binhs, values, errors


def histo_print_per_state(list_of_dict, propname, title, bin=1, bin_histo=10, list_state=None, nadiab=2):
    n = int(len(list_of_dict) / bin)
    llod = [list_of_dict[i:i + n] for i in xrange(0, len(list_of_dict), n)]
    llos = [list_state[i:i + n] for i in xrange(0, len(list_state), n)]  #

    states_dict = {i: [[] for x in xrange(bin)] for i in range(1, nadiab + 1)}

    for i in range(bin):
        for dict, states in zip(llod[i], llos[i]):
            for time in dict:
                if states.get(time) is not None:
                    states_dict[states.get(time)[0]][i] += dict.get(time)

    for state in range(1, nadiab + 1):
        filename = '%s-histo-state-%d-%s.dat' % (propname, state, title)
        file = open(filename, 'w')
        file.write('# %s(%s) Density_of_probability  Error\n' % (propname, units.get(propname)))
        slist = states_dict[state]
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
                bins.append(np.histogram(slist[i], bins=bin_histo, range=(min_, max_), normed=True)[1])
                histos.append(np.histogram(slist[i], bins=bin_histo, range=(min_, max_), normed=True)[0])
        for binh, value, error in zip(np.mean(bins, axis=0),
                                      np.append(np.mean(histos, axis=0), [0]),
                                      np.append(np.std(histos, axis=0, ddof=1), [0])):
            line = '%20.10f %20.10f %20.10f ' % (binh, value, error)
            file.write(line + '\n')
        file.close()
        os.system('mv %s data-%s' % (filename, title))


def rmse(array1, array2):
    rms = 0
    for (x, y) in zip(array1, array2):
        rms += (x - y) ** 2
    return np.sqrt(rms / len(array1))


def mean_over_run(list_of_dict, bin=1):
    n = int(len(list_of_dict) / bin)
    llod = [list_of_dict[i:i + n] for i in xrange(0, len(list_of_dict), n)]
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
    energy_barrier = np.square((reorganization + free_energy)) / (4 * reorganization)
    return energy_barrier


def calculate_marcus_na_rate(coupling, reorganization, free_energy, temperature):
    # Convert everything in atomic units
    coupling = coupling / 27.211399
    reorganization = reorganization / 27.211399
    free_energy = free_energy / 27.211399
    temperature = temperature / 315777.09

    energy_barrier = calculate_barrier(reorganization, free_energy)

    rate = (2 * np.pi) * np.square(coupling) * np.exp(- energy_barrier / temperature) / np.sqrt(
        4 * np.pi * reorganization * temperature)
    rate = rate / 0.024188
    return rate


def rabi_oscillation(coupling, delta_e, delta_t, steps):
    atomic_unit_time = 2.418884326505E-17 / 1.0E-15
    delta_t = delta_t / atomic_unit_time
    amplitude = (4 * coupling ** 2) / (delta_e ** 2 + 4 * coupling ** 2)
    frequence = np.sqrt(delta_e ** 2 + 4 * coupling ** 2)

    pop = [amplitude * (np.sin(frequence * step * delta_t / 2)) ** 2 for step in range(steps + 1)]
    return pop


def expo_constraint(x, b):
    return 0.5 * np.exp(- b * x) + 0.5


def expo_free(x, b, a, c):
    return a * np.exp(- b * x) + c


def log_constraint(x, b):
    return -b * x


def log_free(x, b, a, c):
    return np.log(a * np.exp(- b * x) + c)


def create_file(property, title, text,  label='None', tuple = 'None'):
    if tuple is not 'None':
        filename = '%s-%s.dat' % (property, '-'.join(tuple))
    else:
        filename = '%s-%s.dat' % (property, title)
    file = open(filename, 'w')
    unit = units.get(property)
    if (label == 'Mean'):
        replace = (unit,) * 4
        header = '# Run  Mean (%s)  Std (%s)  Drift (%s/fs)   QMean (%s) \n' % replace
    elif (label == 'Spec'):
        header = headers.get(property)
    elif (label == 'Initial'):
        header = '# Run   Initial value (%s)\n' % unit
    else:
        header = None
    if header is not None:
        file.write(header)
    file.write(text)
    file.close()
    return  filename


def print_dat(array, label, title):
    new = map(list, zip(*array))
    filename = label + '.dat'
    file = open(filename, 'w')
    for line in new:
        file.write('    '.join(map(str, line)) + '\n')
    file.close()
    os.system('mv %s data-%s' % (filename, title))



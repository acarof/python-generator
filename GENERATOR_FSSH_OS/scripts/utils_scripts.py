import re
import numpy as np


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
        self.init = False
        self.hamiltonians = self._read_xyz_file("run-hamilt-1.xyz")
        self.coefficients = self._read_xyz_file("run-coeff-1.xyz")


    def extract(self, property, init=False, **kwargs):
        if init:
            self.init = True
            self.timestep  = float(self.get_input_key(['\TIMESTEP'])[0])
        else:
            self.init = False
        if (property == 'Atoms Number'):
            return self._get_atoms_number()
        elif (property == 'Off-diagonals'):
            return self._extract_off_diagonal_elements()
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
            populations = self._extract_population()
            return self._extract_msd(populations)
        elif property == 'IPR':
            populations = self._extract_population()
            return self._extract_ipr(populations)
        elif property == 'Projected-populations':
            return self._extract_projected_populations()
        elif property == 'Projected-MSD':
            populations = self._extract_projected_populations()
            return self._extract_msd(populations)
        elif property == 'Projected-IPR':
            populations = self._extract_projected_populations()
            return self._extract_ipr(populations)
        elif property == 'MQC-populations':
            return self._extract_mqc_populations()
        elif property == 'MQC-MSD':
            populations = self._extract_mqc_populations()
            return self._extract_msd(populations)
        elif property == 'MQC-IPR':
            populations = self._extract_mqc_populations()
            return self._extract_ipr(populations)
        else:
            print "Extraction of %s not implemented" % property
            raise SystemExit

    def _extract_mqc_populations(self):
        result = {}
        states = self._extract_state()
        for time in self.hamiltonians:
            if (self.coefficients.get(time) is not None) and (states.get(time) is not None):
                diabat = [np.complex(coeff[2], coeff[3]) for coeff in self.coefficients[time]]
                pop = np.array([np.absolute(x) ** 2 for x in diabat])
                #print "populations", pop

                hamilt = np.transpose( self.hamiltonians[time])[2:]
                eigenvalues, eigenvectors = np.linalg.eig( hamilt )
                idx = eigenvalues.argsort()
                eigenvectors = eigenvectors[:, idx]

                state = states[time][0] - 1
                pop += np.array([ np.absolute(x)**2 for x in eigenvectors[:, state]])
                #print "projected", np.array([ np.absolute(x)**2 for x in eigenvectors[:, state]])
                #print "pop", pop

                adiabat = np.dot( eigenvectors.transpose(), diabat)
                pop_adiab = [np.absolute(x) ** 2 for x in adiabat]
                #print "pop_adiab", pop_adiab
                for state in range(len(pop_adiab)):
                    pop += - pop_adiab[state]*np.array([ np.absolute(x)**2 for x in eigenvectors[:, state]])
                    #print "adiab_pop",  pop_adiab[state]
                    #print "eigen", np.array([ np.absolute(x)**2 for x in eigenvectors[:, state]])
                    #print "correction", pop

                result[time] = pop
        return result

    def _extract_projected_populations(self):
        result = {}
        states = self._extract_state()
        for time in self.hamiltonians:
            if states.get(time) is not None:
                hamilt = np.transpose( self.hamiltonians[time])[2:]
                #print hamilt
                eigenvalues, eigenvectors = np.linalg.eig( hamilt )
                #print eigenvalues
                #print eigenvectors
                idx = eigenvalues.argsort()
                eigenvalues = eigenvalues[idx]
                #print eigenvalues
                eigenvectors = eigenvectors[:, idx]
                #print eigenvectors

                state = states[time][0] - 1 # warning
                #print state
                #print eigenvalues
                coeff = eigenvectors[:, state]
                #print coeff
                #print np.dot(hamilt, coeff)
                #print coeff
                populations =  [ np.absolute(x)**2 for x in coeff]
                #print populations
                result[time] = populations
        return result


    def _extract_ipr(self, populations):
        #populations = self._extract_population()
        result = {}
        for time, pop in populations.items():
            ipr = 0
            for popi in pop:
                ipr += popi**2
            ipr = 1.0/ipr
            result[time] = ipr
        return result

    def _extract_msd(self, populations):
        #populations = self._extract_population()
        com = self._extract_com()
        result = {}
        pos0 = 0.0
        for (x,y) in zip(populations[0.0], com[0.0]):
            pos0 +=  x *y
        print "pos0", pos0
        for time, pop in populations.items():
            if com.get(time):
                msd = 0.0
                for (x,y) in zip(pop, com[time]):
                    msd += x*(y)
                msd -= pos0
                result[time] = msd*msd
        return result

    def _extract_com(self, filename = 'run-pos-1.xyz'):
        populations = self._extract_population()
        results = {}
        for time, pop in populations.items():
            results[time] = [ x*6.038 for x in range(len(pop))]
        return results


    def _extract_state(self, filename='run-sh-1.log'):
        try:
            results = self._states
        except:
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
            self._states = results
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
        result = {}
        for time in self.hamiltonians:
            if self.coefficients.get(time) is not None:
                hamilt = np.transpose( self.hamiltonians[time])[2:]
                eigenvalues, eigenvectors = np.linalg.eig( hamilt )
                idx = eigenvalues.argsort()
                eigenvalues = eigenvalues[idx]
                eigenvectors = eigenvectors[:, idx]
                diabat = [ np.complex(coeff[2], coeff[3]) for coeff in self.coefficients[time] ]
                adiabat = np.dot( eigenvectors.transpose(), diabat)
                #adiabat = eigenvectors.transpose().dot(diabat)
                #print "Test", adiabat - adiabat1
                populations = [ np.absolute(x)**2 for x in adiabat]
                result[time] = populations
        return result

    def _extract_adiabatic_energies(self):
        result = {}
        for time in self.hamiltonians:
            hamilt = np.transpose( self.hamiltonians[time])[2:]
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
        delta_e = {}
        for time in self.hamiltonians:
            delta_e[time] = []
            for i in range(1, len(self.hamiltonians.get(time))):
                delta_e[time].append(
                    self.hamiltonians.get(time)[i][i + 2] - self.hamiltonians.get(time)[i - 1][i - 1 + 2]
                )
        return delta_e

    def _extract_coupling(self):
        couplings = {}
        for time in self.hamiltonians:
            couplings[time] = [self.hamiltonians.get(time)[0][3]]
        return couplings

    def _extract_off_diagonal_elements(self):
        nadiab = int(self.get_input_key(['NUMBER_DIABATIC_STATES'])[0])
        couplings = {}
        for time in self.hamiltonians:
            couplings[time] = []
            for state in range(nadiab-1):
                couplings[time].append( self.hamiltonians[time][state][state + 3] )
        return couplings


    def _extract_population(self, filename='run-coeff-1.xyz'):
        populations = {}
        for time in self.coefficients:
            list = self.coefficients.get(time)
            pop = []
            for i, coeff in enumerate(list):
                pop.append(np.square(np.absolute(np.complex(coeff[2], coeff[3]))))
            populations.update({time: pop})
        return populations





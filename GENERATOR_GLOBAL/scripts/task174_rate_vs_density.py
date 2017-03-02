import string, re, struct, sys, math, os, time
import numpy as np
import importlib, imp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from operator import itemgetter
from utils_scripts import *



densities = [0.01, 0.02, 0.03, 0.04]
detail_properties = ['Populations']
histo_properties = ['Delta_E']
mean_properties = ['Couplings', 'Temperature']
specific_properties = ['FSSH']
total_properties = detail_properties + mean_properties + specific_properties + histo_properties

reorganization = 0.300
free_energy = 0.00

density_dict = {
    218 : 0.01,
    506 : 0.02,
    712 : 0.03,
    979 : 0.04
}

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


if not os.path.isdir('data'):
    os.mkdir('data')

name_bucket = os.getcwd().split('/')[-1]


nadiab = 2
os.system('cd ..')
for i, directory in enumerate(os.listdir('.')):
    if i % 100 == 0:
        print time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    if 'run' in directory:
        print "Do %s" % directory
        os.chdir(directory)
        dir = FSSHRun(directory)
        natoms = dir.extract('Atoms Number')
        density = density_dict.get(natoms)
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
                histo = histogram(prop)
                properties_.get(property).dict.get(density).append(histo)
                dir.detailed_print(prop, property)
                dir.histo_print(histo, property)
        os.system('mv *.dat ../data/')
        os.chdir('..')

for property in (mean_properties):
    create_file(property, name_bucket, properties_.get(property).info, 'Mean')
for property in (specific_properties):
    create_file(property, name_bucket, properties_.get(property).info, 'Spec')


plt.figure(1, figsize=(6,5))
rates=[]
marcus_rates = []
final = ''
for density in densities:
    pop = mean_over_run( properties_.get('Populations').dict.get(density) )

    for adiab in range(1, nadiab):
        ax =plt.gca()
        previous, = ax.plot(pop[0], pop[adiab], label = r'%s (atoms/$\AA^3$)' % density)
        try:
            popt, pcov = curve_fit(expo, pop[0], pop[adiab])
        except:
            popt = 0.0, 0.0, 0.0
        rate = popt[1] / 2
        rates.append( rate)
        plt.plot( pop[0], [ expo( x, *popt) for x in pop[0] ], color = previous.get_color(), linestyle = '--' )
        print "The decay rate for state %d is %f fs-1" %(adiab, rate)

    couplings = np.array( properties_.get('Couplings').dict.get(density)).transpose()[3] * 27.211399 * 1000
    print "Couplings values (meV): "
    print couplings
    mean_coupling =  np.mean( couplings )
    print "The average coupling is %f meV " %  mean_coupling

    temperatures = np.array( properties_.get('Temperature').dict.get(density)).transpose()[0]
    print "Temperature values (K): "
    print temperatures
    mean_temperature =  np.mean( temperatures )
    print "The average temperature is %f K " %  mean_temperature

    marcus_rate = calculate_marcus_na_rate(mean_coupling / 1000, reorganization, free_energy, mean_temperature)
    marcus_rates.append( marcus_rate)
    print "The Marcus rate is  %f fs-1" % marcus_rate

    final += '%20.10f  %20.10f %20.10f %20.10f %20.10f \n' % \
             (density, rate, marcus_rate, mean_coupling, mean_temperature)


create_file('Rates', name_bucket, final)
final = 'Density (atom/A^3    Num.Rate (fs-1)    MarcusRate (fs-1)    Coupling (meV)    Temperature (K)\n'
create_file('Rates-info', name_bucket, final)


plt.title('Population evolution for different solvent density')
plt.xlabel('Time (fs)')
plt.locator_params(axis='x', nbins=5)
plt.ylabel('Population')
plt.legend()
plt.legend(bbox_to_anchor=(1.6, 1))
plt.savefig('population_vs_time_%s.png' % name_bucket, bbox_inches='tight')
#plt.ylim([5, 36])
plt.legend(bbox_to_anchor=(2.01, 1))



plt.figure(2, figsize=(6,5))
plt.title('FSSH and Marcus rates vs density')
plt.plot(densities, rates, label = 'FSSH rates', marker = 'o')
plt.plot(densities, marcus_rates, label = 'Marcus rate', marker = 'v')
plt.xlabel(r'Density (atoms/$\AA^3$)')
plt.yscale('log')
plt.ylabel(r'Rate (fs$^{-1}$)')
plt.legend()
plt.legend(bbox_to_anchor=(1.6, 1))
plt.savefig('rate_vs_density_%s.png' % name_bucket, bbox_inches='tight')


os.system('mv *.png data/')

#plt.show()
plt.close()



print "End of the analysis"
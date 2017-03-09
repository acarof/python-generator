import string, re, struct, sys, math, os, time
import numpy as np
import importlib, imp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from operator import itemgetter
from utils_scripts import *



#densities = [0.01, 0.02, 0.03, 0.04]
detail_properties = ['Populations', 'State']
histo_properties = ['Delta_E']
mean_properties = ['Couplings', 'Temperature']
specific_properties = ['FSSH']
total_properties = detail_properties + mean_properties + specific_properties + histo_properties
test = True

if test:
    dirlist = ['run-%d' % i for i in range(40)]
else:
    dirlist = os.listdir('.')

bin = 5
bin_histo = 10

reorganization = 0.300
free_energy = 0.00

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
density_dict =  {v: k for k, v in dict_natoms.iteritems()}
#density_dict = {
#     12 : 0.00,
#    218 : 0.01,
#    506 : 0.02,
#    712 : 0.03,
#    979 : 0.04
#}
dict_sizebox = {
    0.00: [60.0, 60.0, 60.0],
    0.0001: [60.0, 60.0, 60.0],
    0.0005: [60.0, 60.0, 60.0],
    0.001: [60.0, 60.0, 60.0],
    0.005: [60.0, 60.0, 60.0],
    0.01: [30.0, 30.0, 30.0],
    0.02: [30.0, 30.0, 30.0],
    0.03: [30.0, 30.0, 30.0],
    0.04: [30.0, 30.0, 30.0],
    0.06: [20.0, 20.0, 20.0],
    0.08: [20.0, 20.0, 20.0],
    0.10: [20.0, 20.0, 20.0]
}

densities = density_dict.values()
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
        density = density_dict.get(natoms)
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
                #histo = histogram(prop)
                properties_.get(property).dict.get(density).append(prop)
                dir.detailed_print(prop, property)
                #dir.histo_print(histo, property)
        os.system('mv *run*.dat ../per-run/')
        os.chdir('..')


for property in (mean_properties):
    create_file(property, title, properties_.get(property).info, 'Mean')
for property in (specific_properties):
    create_file(property, title, properties_.get(property).info, 'Spec')


plt.figure(1, figsize=(6,5))
mean_rates=[]
mean_rates_c = []
mean_marcus_rates = []
error_rates=[]
error_rates_c = []
error_marcus_rates = []
final = ''
final_c = ''
final_marcus = ''
my_densities.sort()
for density in my_densities:
    natoms = dict_natoms.get(density)
    sizebox = dict_sizebox.get(density)
    volume = np.prod( sizebox)
    real_density = natoms/volume
    print "Start density: %f, real density %f" % (density, real_density)
    for property in histo_properties:
        histo_print(properties_.get(property).dict.get(density), property, title, bin, bin_histo, properties_.get('State').dict.get(density))

    pops = mean_over_run( properties_.get('Populations').dict.get(density), bin  )
    rates = []
    rates_c = []
    marcus_rates = []
    A = []
    B = []
    for i in range(bin):
        print "Start the bin %i" % i
        pop = pops[i]
        print_dat(pop, 'Populations-bin-%d-%s' % (i, title), title )

        for adiab in range(1, nadiab):
            #print "LogPop"
            #for x in pop[adiab]:
            #    print  x, np.log(x - 0.5) - np.log(0.5)
            try:
                popt, pcov = curve_fit(expo_free, pop[0], pop[adiab])
                rate = (popt[0] / 2) * 1E15
                rates.append(rate)
                A.append(popt[1])
                B.append(popt[2])
            except:
                print "Can't fit without constraint for bin %d" % i
                pass
            try:
            #popt_c, pcov_c = curve_fit( log_constraint, pop[0],  [np.log(x - 0.5) - np.log(0.5) for x in pop[adiab]] )
                popt_c, pcov_c = curve_fit( expo_constraint, pop[0], pop[adiab])
                rate_c = popt_c[0] / 2 * 1E15
                rates_c.append(rate_c)
            #popt_c = np.polyfit( pop[0], [np.log(x - 0.5) - np.log(0.5) for x in pop[adiab]], 1 )
            except:
                print "Can't fit with constraint for bin %d" % i
                pass
            #print popt_c
            #sys.exit()

    print "A"
    print A
    print "B"
    print B

    rate = np.mean(rates)
    mean_rates.append(rate)
    error_rates.append(np.std(rates, ddof=1))
    final += '%.3e  %.3e %.3e  %f  %f \n' % \
             (real_density, rate, np.std(rates, ddof=1),  np.mean(A), np.mean(B))

    rate_c = np.mean(rates_c)
    mean_rates_c.append(rate_c)
    error_rates_c.append(np.std(rates_c, ddof=1))
    final_c += '%.3e  %.3e %.3e  \n' % \
             (real_density, rate_c, np.std(rates_c, ddof=1))

    ax = plt.gca()
    ave_pop = np.mean(pops, axis = 0)
    previous, = ax.plot( ave_pop[0], ave_pop[1], label=r'%s (atoms/$\AA^3$)' % density)
    #print rate, np.mean(A), np.mean(B)
    #print [expo_free(x, rate* 2 / 1E15, np.mean(A), np.mean(B) ) for x in ave_pop[0] ]
    #sys.exit()
    plt.plot(ave_pop[0], [expo_free(x, rate* 2 / 1E15, np.mean(A), np.mean(B) ) for x in ave_pop[0] ], color=previous.get_color(), linestyle='--')
    plt.plot(ave_pop[0], [expo_constraint(x, rate_c* 2 / 1E15) for x in ave_pop[0]], color=previous.get_color(), linestyle=':')


    couplings = np.array( properties_.get('Couplings').dict.get(density)).transpose()[3] * 27.211399 * 1000
    temperatures = np.array(properties_.get('Temperature').dict.get(density)).transpose()[0]
    n = int(len(couplings) / bin)
    lcouplings = [couplings[i:i + n] for i in xrange(0, len(couplings), n)]
    ltemperatures = [temperatures[i:i  + n] for i in xrange(0, len(temperatures), n)]
    std_temperatures = []
    for i in range(bin):
        mean_coupling = np.mean(lcouplings[i])
        #print "Couplings values (meV): ", lcouplings[i]
        #print "The average coupling is %f meV " % mean_coupling

        mean_temperature = np.mean(ltemperatures[i])
        std_temperatures.append( np.std(ltemperatures[i]) )
        #print "Temperature values (K): ", ltemperatures[i]
        #print "The average temperature is %f K " % mean_temperature

        marcus_rate = calculate_marcus_na_rate(mean_coupling / 1000, reorganization, free_energy, mean_temperature) * 1E15
        marcus_rates.append( marcus_rate)
        #print "The Marcus rate is  %f fs-1" % marcus_rate

    marcus_rate = np.mean( marcus_rates)
    mean_marcus_rates.append(marcus_rate)
    error_marcus_rates.append( np.std(marcus_rates, ddof=1))

    final_marcus += '%.3e  %.3e %.3e  %f %f %f\n' % \
             (real_density,
              marcus_rate, np.std(marcus_rates, ddof=1), mean_coupling, mean_temperature, np.mean(std_temperatures))


create_file('Rates-without-constraint', title, final)
create_file('Rates-with-constraint', title, final_c)
create_file('Rates-Marcus', title, final_marcus)
final = 'Density (atom/A^3    Num.Rate (s-1) Error  A  B  \n'
final_c = 'Density (atom/A^3      Num.RateConstraint (s-1)  Error MarcusRate (s-1) Error \n'
final_marcus = 'Density (atom/A^3    MarcusRate (s-1) Error   Coupling (meV)    Temperature (K)\n'
create_file('Rates-info-without-constraint', title, final)
create_file('Rates-info-with-constraint', title, final_c)
create_file('Rates-info-Marcus', title, final_marcus)

short_time = time.strftime("%y%m%d", time.localtime())
plt.title('Population evolution for different solvent density')
plt.xlabel('Time (fs)')
plt.locator_params(axis='x', nbins=5)
plt.ylabel('Population')
plt.legend()
plt.legend(bbox_to_anchor=(1.6, 1))
plt.savefig('population_vs_time_%s.png' % title, bbox_inches='tight')
#plt.ylim([5, 36])
plt.legend(bbox_to_anchor=(2.01, 1))


plt.figure(2, figsize=(6,5))
plt.title('FSSH and Marcus rates vs density')
plt.errorbar(my_densities, mean_rates, yerr = error_rates, label = 'FSSH rates', marker = 'o')
plt.errorbar(my_densities, mean_rates_c, yerr = error_rates_c, label = 'FSSH rates (cons)', marker = 'x')
plt.errorbar(my_densities, mean_marcus_rates, yerr = error_marcus_rates, label = 'Marcus rate', marker = 'v')
plt.xlabel(r'Density (atoms/$\AA^3$)')
plt.yscale('log')
plt.ylabel(r'Rate (fs$^{-1}$)')
plt.legend()
plt.legend(bbox_to_anchor=(1.6, 1))
plt.savefig('rate_vs_density_%s.png' % title, bbox_inches='tight')


os.system('mv *.png %s' % dataname)

#plt.show()
plt.close()



print "End of the analysis"
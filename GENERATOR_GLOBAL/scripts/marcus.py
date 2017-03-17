import string, re, struct, sys, math, os, time
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad




def landau_energy_diabatic_a( delta_e, reorga, free_energy ):
    return ( 1 / (4* reorga)) * ( delta_e - ( free_energy + reorga))**2



def landau_energy_diabatic_b( delta_e, reorga, free_energy ):
    return ( 1 / (4* reorga)) * ( delta_e - ( free_energy - reorga))**2 + free_energy



def landau_energy_adiabatic_ground( delta_e, reorga, free_energy, coupling):
    first = (landau_energy_diabatic_a(delta_e, reorga, free_energy) + landau_energy_diabatic_b(delta_e, reorga, free_energy))/ 2.0
    second = np.sqrt( delta_e**2 + 4 * coupling**2) / 2.0
    return first - second



def landau_energy_adiabatic_excited( delta_e, reorga, free_energy, coupling):
    first = (landau_energy_diabatic_a(delta_e, reorga, free_energy) + landau_energy_diabatic_b(delta_e, reorga, free_energy))/ 2.0
    second = np.sqrt( delta_e**2 + 4 * coupling**2) / 2.0
    return first + second



def boltzmann_weight( delta_e, landau_energy,  reorga, free_energy, coupling, temperature ):
    k_b = 8.6173303 * 10**(-5) # Boltzmann constant in eV/K
    beta = 1 / (k_b * temperature)
    return np.exp( - beta * landau_energy(delta_e, reorga, free_energy, coupling))



def average_dict(dict1, number):
    result = {}
    for key in dict1:
        result[key] = np.array( dict1[key] ) / number
    return result



def calculate_integral_boltzmann_weight(landau_energy, reorga, free_energy, coupling, temperature):
    return quad( boltzmann_weight, - np.inf, np.inf, args = (landau_energy, reorga, free_energy, coupling, temperature)  )[0]



def calculate_boltzman_ratio( reorga, free_energy, coupling, temperature ):
    A0 = calculate_integral_boltzmann_weight(landau_energy_adiabatic_ground, reorga, free_energy, coupling, temperature)
    A1 = calculate_integral_boltzmann_weight(landau_energy_adiabatic_excited, reorga, free_energy, coupling, temperature)
    return A0 / ( A1 + A0)






#reorga = 1.0
#free_energy = -0.1
#coupling = 0.07
#temperature = 300




#print calculate_integral_boltzmann_weight(landau_energy_adiabatic_ground, reorga, free_energy, coupling, temperature)
#
#print calculate_boltzman_ratio(reorga, free_energy, coupling, temperature )
#
#
#
#delta_e_axis = np.arange( -2, 2, 0.1)
##sys.exit()
#file = open('test.dat', 'w')
#for delta_e in delta_e_axis:
#    file.write( '%f    %f\n' % (delta_e, boltzmann_weight( delta_e, landau_energy_adiabatic_ground,  reorga, free_energy, coupling, temperature ))   )
#file.close()
#file = open('test2.dat', 'w')
#for delta_e in delta_e_axis:
#    file.write( '%f    %f\n' % (delta_e, boltzmann_weight( delta_e, landau_energy_adiabatic_excited,  reorga, free_energy, coupling, temperature ))   )
#file.close()
#
#
#plt.plot( delta_e_axis, [ boltzmann_weight( delta_e, landau_energy_adiabatic_ground, reorga, free_energy, coupling, temperature ) for delta_e in delta_e_axis ])
#plt.plot( delta_e_axis, [ boltzmann_weight( delta_e, landau_energy_adiabatic_excited,  reorga, free_energy, coupling, temperature ) for delta_e in delta_e_axis ])
#plt.yscale('log')
#plt.show()
#sys.exit()
#
#
#plt.plot( delta_e_axis, [ landau_energy_diabatic_a( delta_e, reorga, free_energy ) for delta_e in delta_e_axis ])
#plt.plot( delta_e_axis, [ landau_energy_diabatic_b( delta_e, reorga, free_energy ) for delta_e in delta_e_axis ])
#plt.plot( delta_e_axis, [ landau_energy_adiabatic_ground( delta_e, reorga, free_energy, coupling ) for delta_e in delta_e_axis ])
#plt.plot( delta_e_axis, [ landau_energy_adiabatic_excited( delta_e, reorga, free_energy, coupling ) for delta_e in delta_e_axis ])
#plt.show()
#
#sys.exit()
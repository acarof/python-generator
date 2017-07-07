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



def boltzmann_weight( delta_e, landau_energy, temperature, *args ):
    k_b = 8.6173303 * 10**(-5) # Boltzmann constant in eV/K
    beta = 1 / (k_b * temperature)
    return np.exp( - beta * landau_energy(delta_e, *args))


def calculate_integral_boltzmann_weight(landau_energy, temperature, *args):
    return quad( boltzmann_weight, - np.inf, np.inf, args = (landau_energy, temperature,)  + args )[0]

def calculate_free_energy(landau_energy, temperature, *args ):
    I = calculate_integral_boltzmann_weight(landau_energy, temperature, *args)
    k_b = 8.6173303 * 10**(-5) # Boltzmann constant in eV/K
    return - k_b * temperature * np.log( I )



def calculate_boltzman_ratio( reorga, free_energy, coupling, temperature, state ):
    I0 = calculate_integral_boltzmann_weight(landau_energy_adiabatic_ground, temperature, reorga, free_energy, coupling, )
    I1 = calculate_integral_boltzmann_weight(landau_energy_adiabatic_excited, temperature, reorga, free_energy, coupling)
    if state == 'Excited':
        return I1 / (I0 + I1)
    else:
        return I0 / ( I1 + I0)




def calculate_barrier(reorganization, free_energy):
    energy_barrier = np.square( ( reorganization + free_energy ) ) / ( 4 * reorganization)
    print energy_barrier
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


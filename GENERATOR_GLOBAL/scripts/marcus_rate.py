import string, re, struct, sys, math, os, time
import numpy as np
#import importlib, imp
#import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit

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


coupling = float(raw_input("Electronic couplings (eV): "))
reorganization = float(raw_input("Reorganization energy (eV): ") )
free_energy = float(raw_input("Free energy differences (eV): ") )
temperature = float(raw_input("Temperature (K): ") )

rate = calculate_marcus_na_rate(coupling, reorganization, free_energy, temperature)

print "The rate is %s fs^-1 " % rate


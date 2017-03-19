from marcus import *
from utils_scripts import *
import matplotlib.pyplot as plt





reorga = 1.0
free_energy = -0.1
coupling = 0.07
temperature = 300
k_b =  8.6173303 * 10**(-5) # Boltzmann constant in eV/K
beta = 1 / (k_b * temperature)
delta_e_axis = np.arange( -2, 2, 0.1)
hab_axis = np.arange(0, 0.300, 0.0001)

#plt.show()

# CHECK LANDAU ENERGIES - COMPARISON WITH JOCHEN'S GRAPH

plt.plot( delta_e_axis, [ landau_energy_diabatic_a( delta_e, reorga, free_energy ) for delta_e in delta_e_axis ])
plt.plot( delta_e_axis, [ landau_energy_diabatic_b( delta_e, reorga, free_energy ) for delta_e in delta_e_axis ])
plt.plot( delta_e_axis, [ landau_energy_adiabatic_ground( delta_e, reorga, free_energy, coupling ) for delta_e in delta_e_axis ])
plt.plot( delta_e_axis, [ landau_energy_adiabatic_excited( delta_e, reorga, free_energy, coupling ) for delta_e in delta_e_axis ])
plt.ylabel('Free energy difference (eV)')
plt.xlabel('Vertical energy gap (eV)')
plt.savefig('landau_curves_for_marcus.png')

sys.exit()

# CHECK INTEGRAL: print weight and use xmgrace to check integral value


#sys.exit()
file = open('test.dat', 'w')
for delta_e in delta_e_axis:
    file.write( '%f    %f\n' % (delta_e, boltzmann_weight( delta_e, landau_energy_adiabatic_ground,  reorga, free_energy, coupling, temperature ))   )
file.close()
file = open('test2.dat', 'w')
for delta_e in delta_e_axis:
    file.write( '%f    %f\n' % (delta_e, boltzmann_weight( delta_e, landau_energy_adiabatic_excited,  reorga, free_energy, coupling, temperature ))   )
file.close()



# CHECK INTEGRAL: calculate integral and compare with analytic prediction for a harmonic landau energy (diabatic)
print calculate_integral_boltzmann_weight(landau_energy_diabatic_a, temperature, reorga, free_energy)

print calculate_free_energy(landau_energy_diabatic_a, temperature, reorga, free_energy)




free_energy_diff = np.array([
    calculate_free_energy(landau_energy_adiabatic_excited, temperature, reorga, free_energy, hab) -
    calculate_free_energy(landau_energy_adiabatic_ground, temperature, reorga, free_energy, hab)
    for hab in hab_axis
])

plt.figure(1)
plt.plot( hab_axis, free_energy_diff )
plt.title('Adiabatic free energy difference vs coupling')
plt.xlabel('Coupling (eV)')
plt.ylabel('Free energy difference (eV)')
plt.savefig('adiabat_free_diff_vs_coupling.png')

plt.figure(2)
plt.plot( hab_axis, np.exp( - beta * free_energy_diff) )
plt.title('Ratio adiabatic population vs coupling')
plt.xlabel('Coupling (eV)')
plt.ylabel('Ratio')
plt.savefig('ratio_vs_coupling.png')
import matplotlib.pyplot as plt
from plots import plot_from_file, plot_map

print('rysowanie wykresów ...')

plt.figure(figsize=(6, 6))
plot_from_file('files/glob_1.txt', 'ω = 0.6')
plot_from_file('files/glob_2.txt', 'ω = 1.0')

plt.xlabel('Numer Iteracji')
plt.ylabel('S')
plt.xscale('log')
plt.legend()
plt.suptitle('Relaksacja globalna, S(it)')
plt.show()

omegas = [0.6, 1.0]

file_paths_potential = ['files/map_omega_1.txt', 'files/map_omega_2.txt']
plot_map(file_paths_potential, 'Mapa zrelaksowanego potencjału, ', omegas)

file_paths_err = ['files/err_1.txt', 'files/err_2.txt']
plot_map(file_paths_err, 'Mapa błędu rozwiązania, ', omegas)

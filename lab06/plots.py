import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

def plot_map(file_path, s, var, e2):
    data = np.loadtxt(file_path)
    x_values, y_values, V_values = data[:, 0], data[:, 1], data[:, 2]

    plt.figure(figsize=(7, 6))

    if s == 1:
        max = 10
        min = -10
    else:
        max = 0.8
        min = -0.8

    scatter = plt.scatter(x_values, y_values, c=V_values, cmap=plt.get_cmap('RdBu_r'), vmin = min, vmax = max, marker='s', label='')

    plt.colorbar(scatter, label='V(x, y)')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'nx=ny= {var}, ε1= 1, ε2={e2}')

    plt.tight_layout()
    # plt.show()

file_paths_rho_0 = ['files/map_nx_ny_50.txt', 'files/map_nx_ny_100.txt', 'files/map_nx_ny_200.txt']
nx_ny = [50, 100, 200]
file_paths_rho_n_0 = ['files/map_eps_1_eps_1.txt', 'files/map_eps_1_eps_2.txt', 'files/map_eps_1_eps_10.txt']
e2 = [1, 2, 10]

i=0
for el in file_paths_rho_0:
    plot_map(el, 1, nx_ny[i], 1)
    i+=1

i=0
for el in file_paths_rho_n_0:
    plot_map(el, 2, 100, e2[i])
    i+=1

plt.show()

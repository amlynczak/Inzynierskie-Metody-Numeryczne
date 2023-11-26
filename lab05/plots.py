import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.interpolate import griddata

def plot_map(file_path, k):
    data = np.loadtxt(file_path)
    x_values, y_values, V_values = data[:, 0], data[:, 1], data[:, 2]

    #dane do interpolacji wartosci dla wykresow
    x_range = np.linspace(min(x_values), max(x_values), 100)
    y_range = np.linspace(min(y_values), max(y_values), 100)
    X, Y = np.meshgrid(x_range, y_range)
    V_interpolated = griddata((x_values, y_values), V_values, (X, Y), method='nearest')

    plt.figure(figsize=(6, 6))

    #OG dane
    scatter = plt.scatter(x_values, y_values, c=V_values, cmap='plasma', marker='s', label='Dane oryginalne')

    # Wykres interpolacji
    plt.imshow(V_interpolated, extent=(min(x_values), max(x_values), min(y_values), max(y_values)),
               origin='lower', cmap='plasma', alpha=0.5, aspect='auto')

    plt.colorbar(scatter, label='V(x, y)')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'k = {k}')

    plt.ylim(0, 25)
    plt.xlim(0, 25)

    plt.tight_layout()
    plt.legend()
    #plt.show()

print('rysowanie wykresow...')

file_paths = ['files/map_k_16.txt', 'files/map_k_8.txt', 'files/map_k_4.txt', 'files/map_k_2.txt', 'files/map_k_1.txt']
k = [16, 8, 4, 2, 1]

for i in range(len(file_paths)):
    plot_map(file_paths[i], k[i])

data_S = np.loadtxt('files/s_it.txt')
data_it = np.loadtxt('files/it.txt')

k_16 = []
k_8 = []
k_4 = []
k_2 = []
k_1 = []

k_it = [int(data_it[0]), 
        int(data_it[0]) + int(data_it[1]), 
        int(data_it[0]) + int(data_it[1]) + int(data_it[2]), 
        int(data_it[0]) + int(data_it[1]) + int(data_it[2]) + int(data_it[3]), 
        int(data_it[0]) + int(data_it[1]) + int(data_it[2]) + int(data_it[3]) +  int(data_it[4])]

for i in range(0, k_it[0]):
    k_16.append(data_S[i][1])

for i in range(k_it[0]+1, k_it[1]):
    k_8.append(data_S[i][1])

for i in range(k_it[1]+1, k_it[2]):
    k_4.append(data_S[i][1])

for i in range(k_it[2]+1, k_it[3]):
    k_2.append(data_S[i][1])

for i in range(k_it[3]+1, k_it[4]):
    k_1.append(data_S[i][1])

plt.figure(figsize=(6,6))

plt.plot(range(0, k_it[0]), k_16, label='k = 16', color='blue')
plt.plot(range(k_it[0]+1, k_it[1]), k_8, label='k = 8', color='red')
plt.plot(range(k_it[1]+1, k_it[2]), k_4, label='k = 4', color='yellow')
plt.plot(range(k_it[2]+1, k_it[3]), k_2, label='k = 2', color='black')
plt.plot(range(k_it[3]+1, k_it[4]), k_1, label='k = 1', color='green')
plt.legend(loc='upper right')
plt.xlabel('iteracja')
plt.ylabel('S')

plt.tight_layout()
plt.show()

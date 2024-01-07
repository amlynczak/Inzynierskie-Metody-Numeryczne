import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

def plot_map(file_path, alpha, beta):
    data = np.loadtxt(file_path)
    t_values, x_values, u_values = data[:, 0], data[:, 1], data[:, 2]
    title = f"α = {alpha}, β = {beta}"

    plt.figure(figsize=(7, 3))

    scatter = plt.scatter(t_values, x_values, c=u_values, cmap=plt.get_cmap('plasma'), marker='s', label='')

    plt.colorbar(scatter, label='u(x, t)')

    plt.xlabel('t')
    plt.ylabel('x')
    plt.title(title)

    plt.tight_layout()

alpha = [0.0, 0.0, 0.0, 1.0]
beta = [0.0, 0.1, 1.0, 1.0]

file_paths = ['files/map_a0.0_b0.0.txt', 'files/map_a0.0_b0.1.txt', 'files/map_a0.0_b1.0.txt', 'files/map_a1.0_b1.0.txt']

i=0
for el in file_paths:
    plot_map(file_paths[i], alpha[i], beta[i])
    i+=1

plt.show()

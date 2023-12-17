import matplotlib.pyplot as plt
import math
import numpy as np

def plot_map(file_path):
    data = np.loadtxt(file_path)
    x_values, y_values, V_values = data[:, 0], data[:, 1], data[:, 2]

    plt.figure(figsize=(14, 4))

    scatter = plt.scatter(x_values, y_values, c=V_values, cmap='viridis', marker='s', label='')

    plt.colorbar(scatter, label='V(x, y)')

    plt.xlabel('X')
    plt.ylabel('Y')

    plt.tight_layout()
    #plt.show()


plot_map("files/Vx.txt")
plot_map("files/Vy.txt")

plt.show()

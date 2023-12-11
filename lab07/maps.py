import numpy as np
import matplotlib.pyplot as plt

def read_file(filename):
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    u = data[:, 4]
    v = data[:, 5]

    return x, y, u, v

def plot_map(filename, min1, max1, min2, max2):
    x, y, u, v = read_file(filename)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    
    for i in range(2):
        if i==0:
            scatter = axes[i].scatter(x, y, c=u, cmap='plasma', marker = 's',  vmin=min1, vmax=max1)
            plt.colorbar(scatter, ax=axes[i], label='u(x, y)')
            if filename == "files/Q_-1000.txt":
                axes[i].set_title("Q=-1000, u(x,y) - składowa pozioma")
            else:
                axes[i].set_title("Q=-4000, u(x,y) - składowa pozioma")
        elif i==1:
            scatter = axes[i].scatter(x, y, c=v, cmap='plasma', marker = 's',  vmin=min2, vmax=max2)
            plt.colorbar(scatter, ax=axes[i], label='v(x,y)')
            if filename == "files/Q_-1000.txt":
                axes[i].set_title("Q=-1000, v(x,y) - składowa pionowa")
            else:
                axes[i].set_title("Q=-4000, v(x,y) - składowa pionowa")

        axes[i].set_xlabel('X')
        axes[i].set_ylabel('Y')

    plt.tight_layout()
    #plt.show()

    

plot_map("files/Q_-1000.txt", -2, 16, -6, 1)

plot_map("files/Q_-4000.txt", -10, 70, -14, 4)

plt.show()

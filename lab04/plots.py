import matplotlib.pyplot as plt
import math

def plot_from_file(filename, label):
    data = []

    with open(filename, 'r') as file:
        for line in file:
            iteration, value = map(float, line.split())
            data.append((iteration, value))

    iterations, values = zip(*data)
    label_string = label + ', ' + str(math.floor(iterations[-1])) + ' it'
    plt.plot(iterations, values, label=label_string)

def plot_map(file_paths, name, omegas):
    fig, axes = plt.subplots(1, len(file_paths), figsize=(12, 6))

    for i, file_path in enumerate(file_paths):
        data = []
        with open(file_path, 'r') as file:
            for line in file:
                x, y, err = map(float, line.split())
                data.append((x, y, err))

        x_values, y_values, err_values = zip(*data)

        scatter = axes[i].scatter(x_values, y_values, c=err_values, cmap='viridis', marker='o', edgecolors='black')

        if name == 'Mapa zrelaksowanego potencjału, ':
            plt.colorbar(scatter, ax=axes[i], label='V(x, y)')
        elif name == 'Mapa błędu rozwiązania, ':
            plt.colorbar(scatter, ax=axes[i], label='δ = ∇2V (x, y) + ρ(x, y)/ε')

        axes[i].set_xlabel('X')
        axes[i].set_ylabel('Y')
        axes[i].set_title(name + f'ω = {omegas[i]}')

        axes[i].set_ylim(0, 10)
        axes[i].set_xlim(0, 15)

    plt.tight_layout()

    plt.show()

import matplotlib.pyplot as plt
import numpy as np
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

plot_from_file('files/c_0.txt', 'D=0.0')
plot_from_file('files/c_0.1.txt', 'D=0.1')

plt.xlabel('t_n')
plt.ylabel('c')
plt.legend()
plt.suptitle('c(t_n)')
plt.show()
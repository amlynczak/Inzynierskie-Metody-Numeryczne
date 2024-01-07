from cProfile import label
import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt('files/E_a0.0_b0.0.txt')
data2 = np.loadtxt('files/E_a0.0_b0.1.txt')
data3 = np.loadtxt('files/E_a0.0_b1.0.txt')

t_val_1, E_val_1 = data1[:, 0], data1[:, 1]
t_val_2, E_val_2 = data2[:, 0], data2[:, 1]
t_val_3, E_val_3 = data3[:, 0], data3[:, 1]

# punkt 3, wykresy E(t)
plt.figure(figsize=(8, 6))
plt.xlabel('t')
plt.ylabel('E')
plt.title("E(t)")

plt.plot(t_val_1, E_val_1, label='α=0.0, β=0.0', color='tab:red')
plt.plot(t_val_2, E_val_2, label='α=0.0, β=0.1', color='tab:blue')
plt.plot(t_val_3, E_val_3, label='α=0.0, β=1.0', color='tab:green')

plt.legend(loc='center right') 

#4 punkt, wykres E(t)
data = np.loadtxt('files/E_a1.0_b1.0.txt')
t_val, E_val = data[:, 0], data[:, 1]

plt.figure(figsize=(8, 6))
plt.xlabel('t')
plt.ylabel('E')
plt.title("E(t)")

plt.plot(t_val, E_val, label='α=1.0, β=1.0')

plt.legend(loc='lower right')


plt.show()

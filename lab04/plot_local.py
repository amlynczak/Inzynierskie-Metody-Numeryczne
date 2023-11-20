import matplotlib.pyplot as plt
from plots import plot_from_file

print('rysowanie wykresów ...')

plot_from_file('files/lok_1.txt', 'ω = 1.0')
plot_from_file('files/lok_2.txt', 'ω = 1.4')
plot_from_file('files/lok_3.txt', 'ω = 1.8')
plot_from_file('files/lok_4.txt', 'ω = 1.9')

plt.xlabel('Numer Iteracji')
plt.ylabel('S')
plt.xscale('log')
plt.legend()
plt.suptitle('Relaksacja lokalna, S(it)')
plt.show()

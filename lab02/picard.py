import numpy as np
import matplotlib.pyplot as plt

from  parameters import beta, N, gamma, tmax, dt, u0, TOL, max_iterations, f 

# Metoda Picarda
def picard_method(dt, u0, TOL, max_iterations):
    t_values = [0] #wektor dla t, w ktorych obliczane sa wartosci u(t) oraz z(t)
    u_values = [u0] #wektor dla liczby nosicieli choroby, inicjalizowany u0
    z_values = [N - u0] #wektor dla zdrowych, liczony jako pozostali

    for i in range(int(tmax / dt)):
        t = t_values[-1] #aktualizacja parametrow na podstawie ostatnich elementow
        u = u_values[-1]

        i = 0 #do iterowania
        u_new = u
        while i < max_iterations:
            u_temp = u_new
            u_new = u + (dt / 2) * (f(t, u) + f(t + dt, u_new)) #wzor (9)
            i = i + 1
            if abs(u_new - u_temp) < TOL: #warunek stopu
                break

        t_values.append(t + dt)
        u_values.append(u_new)
        z_values.append(N - u_new)

    return t_values, u_values, z_values

t_values, u_values, z_values = picard_method(dt, u0, TOL, max_iterations) #użycie metody

plt.figure(figsize=(10, 6))
plt.plot(t_values, u_values, label='u(t)')
plt.plot(t_values, z_values, label='z(t)')
plt.xlabel('t')
plt.ylabel('u(t), z(t)')
plt.legend()
plt.grid(True)
plt.title('Metoda Picarda')
plt.show()
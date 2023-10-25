import numpy as np
import matplotlib.pyplot as plt

# Dane wejsciowe
t_max = 5
y0 = 1
lambda_val = -1
delta_t_val = [0.01, 0.1, 1.0]

# rownanie, ktore chcemy rozwiązać f(t, y) - w naszym przypadku przyjmuje dodatkowo wartosc lambdy
def f(t, y, lambda_val):
    return lambda_val * y #(1) wzor, f(t, y) = lam * y

### -- Metoda jadna Eulera -- ## przyjuje parametry potrzebne do obliczenia ze wzoru (3), zwraca tablice z punktami oraz tablice z wartosciami w tych punktach
def euler_method(lambda_val, t_max, delta_t):
    num_steps = int(t_max / delta_t)
    t_values = np.linspace(0, t_max, num_steps + 1) 
    y_values = np.zeros(num_steps + 1)
    y_values[0] = y0

    for i in range(num_steps):
        y_values[i + 1] = y_values[i] + delta_t * f(t_values[i], y_values[i], lambda_val) #wzor (3)

    return t_values, y_values

## -- Metoda jawna RK2 (trapezow) -- ## przyjumuje parametry potrzebne do obliczenia ze wzorow (4), (5), oraz (6), zwraca takie same informacje jak euler
def rk2_method(lambda_val, t_max, delta_t):
    num_steps = int(t_max / delta_t)
    t_values = np.linspace(0, t_max, num_steps + 1)
    y_values = np.zeros(num_steps + 1)
    y_values[0] = y0

    #obliczanie wartosci w punktach zgodnie z instrukcja
    for i in range(num_steps):
        k1 = f(t_values[i], y_values[i], lambda_val)
        k2 = f(t_values[i] + delta_t, y_values[i] + delta_t * k1, lambda_val)
        y_values[i + 1] = y_values[i] + (delta_t / 2) * (k1 + k2)

    return t_values, y_values

## -- Metoda jawna RK4 -- ## wykorzystanie wzorow (7) - (11)
def rk4_method(lambda_val, t_max, delta_t):
    num_steps = int(t_max / delta_t)
    t_values = np.linspace(0, t_max, num_steps + 1)
    y_values = np.zeros(num_steps + 1)
    y_values[0] = y0

    #obliczanie wartosci zgodnie z instrukcja
    for i in range(num_steps):
        k1 = f(t_values[i], y_values[i], lambda_val)
        k2 = f(t_values[i] + delta_t / 2, y_values[i] + (delta_t / 2) * k1, lambda_val)
        k3 = f(t_values[i] + delta_t / 2, y_values[i] + (delta_t / 2) * k2, lambda_val)
        k4 = f(t_values[i] + delta_t, y_values[i] + delta_t * k3, lambda_val)
        y_values[i + 1] = y_values[i] + (delta_t / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return t_values, y_values

#funkcja obliczajaca bald globalny - wykorzystuje pozniej przy tworzeniu wykresow
def global_error(analytical, numerical):
    return (numerical - analytical)


plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
for delta_t in delta_t_val:
    t_values, y_values = euler_method(lambda_val, t_max, delta_t)
    analytical_values = np.exp(lambda_val * t_values)
    error_values = global_error(analytical_values, y_values)
    plt.plot(t_values, y_values, 'o', label=f'delta_t = {delta_t}')
plt.plot(t_values, analytical_values, 'k--', label='analityczne')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.legend()
plt.title('z.1 - Metoda Eulera - rozwiązanie')

plt.subplot(1, 2, 2)
for delta_t in delta_t_val:
    t_values, y_values = euler_method(lambda_val, t_max, delta_t)
    analytical_values = np.exp(lambda_val * t_values)
    error_values = global_error(analytical_values, y_values)
    plt.plot(t_values, error_values, label=f'delta_t = {delta_t}')
plt.xlabel('t')
plt.ylabel('y_num(t) - y_dok(t)')
plt.legend()
plt.title('z.1 - Metoda Eulera - błąd globalny')
#plt.show()

# Metoda RK2 (Trapezów)
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
for delta_t in delta_t_val:
    t_values, y_values = rk2_method(lambda_val, t_max, delta_t)
    analytical_values = np.exp(lambda_val * t_values)
    error_values = global_error(analytical_values, y_values)
    plt.plot(t_values, y_values, 'o', label=f'delta_t = {delta_t}')
plt.plot(t_values, analytical_values, 'k--', label='analityczne')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.legend()
plt.title('z.2 - Metoda RK2 - rozwiązanie')

plt.subplot(1, 2, 2)
for delta_t in delta_t_val:
    t_values, y_values = rk2_method(lambda_val, t_max, delta_t)
    analytical_values = np.exp(lambda_val * t_values)
    error_values = global_error(analytical_values, y_values)
    plt.plot(t_values, error_values, label=f'delta_t = {delta_t}')
plt.xlabel('t')
plt.ylabel('y_num(t) - y_dok(t)')
plt.legend()
plt.title('z.2 - Metoda RK2 - błąd globalny')
#plt.show()

# Metoda RK4
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
for delta_t in delta_t_val:
    t_values, y_values = rk4_method(lambda_val, t_max, delta_t)
    analytical_values = np.exp(lambda_val * t_values)
    error_values = global_error(analytical_values, y_values)
    plt.plot(t_values, y_values, 'o', label=f'delta_t = {delta_t}')
plt.plot(t_values, analytical_values, 'k--', label='Rozwiązanie analityczne')
plt.xlabel('t')
plt.ylabel('y(t)')
plt.legend()
plt.title('z.3 - Metoda RK4 - rozwiązanie')

plt.subplot(1, 2, 2)
for delta_t in delta_t_val:
    t_values, y_values = rk4_method(lambda_val, t_max, delta_t)
    analytical_values = np.exp(lambda_val * t_values)
    error_values = global_error(analytical_values, y_values)
    plt.plot(t_values, error_values, label=f'delta_t = {delta_t}')
plt.xlabel('t')
plt.ylabel('y_num(t) - y_dok(t)')
plt.legend()
plt.title('z.3 - Metoda RK4 - błąd globalny')
plt.show()
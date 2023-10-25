import numpy as np
import matplotlib.pyplot as plt

# dane wejsciowe dla parametrow
dt = 1e-4
R = 100
L = 0.1
C = 0.001
omega0 = np.sqrt(1 / (L * C))
T0 = 2 * np.pi / omega0
t_end = 4 * T0
t_values = np.arange(0, t_end, dt)
Q0 = 0
I0 = 0

# potencjal zrodla napiecia V(t) (dodatkowo omega_V)
def V(t, omega_V):
    return 10 * np.sin(omega_V * t)

# wzor (13)
def f(t, Q, I):
    return I

#wzor (14)
def g(t, Q, I, omegaV):
    tmp1 = V(t, omegaV)/L
    tmp2 = R/L * I
    tmp3 = 1/(L*C) * Q
    return tmp1 - tmp2 - tmp3

# tablice wynikowe
Q_values = np.zeros(len(t_values))
I_values = np.zeros(len(t_values))

# tablica ze zmiennymi omegaV
omegaV_val = [0.5 * omega0, 0.8 * omega0, omega0, 1.2 * omega0]
for omegaV in omegaV_val:
    Q = Q0
    I = I0

    #obliczenie elementow k.. do koncowych wzorow oraz obliczenie samych wartosci pradu i ladunku
    for i in range(len(t_values)):
        kQ1 = f(t_values[i], Q, I)
        kI1 = g(t_values[i], Q, I, omegaV)

        kQ2 = f(t_values[i] + dt / 2, Q + dt / 2 * kQ1, I + dt / 2 * kI1)
        kI2 = g(t_values[i] + dt / 2, Q + dt / 2 * kQ1, I + dt / 2 * kI1, omegaV)

        kQ3 = f(t_values[i] + dt / 2, Q + dt / 2 * kQ2, I + dt / 2 * kI2)
        kI3 = g(t_values[i] + dt / 2, Q + dt / 2 * kQ2, I + dt / 2 * kI2, omegaV)

        kQ4 = f(t_values[i] + dt, Q + dt * kQ3, I + dt * kI3)
        kI4 = g(t_values[i] + dt, Q + dt * kQ3, I + dt * kI3, omegaV)

        Q = Q + dt / 6 * (kQ1 + 2 * kQ2 + 2 * kQ3 + kQ4)
        I = I + dt / 6 * (kI1 + 2 * kI2 + 2 * kI3 + kI4)

        Q_values[i] = Q
        I_values[i] = I

    # wykres dla Q(t) dla danego omega_V
    plt.figure(1)
    plt.plot(t_values, Q_values, label=f'ωV/ω0 = {omegaV/omega0:.1f}')
    plt.xlabel('t')
    plt.ylabel('Q(t)')
    plt.title('z.4 - Metoda RK4, Q(t)')
    plt.grid(True)
    plt.legend()

    # wykres dla I(t) dla danego omega_V
    plt.figure(2)
    plt.plot(t_values, I_values, label=f'ωV/ω0 = {omegaV/omega0:.1f}')
    plt.xlabel('x')
    plt.ylabel('I(t)')
    plt.title('z.4 - Metoda RK4 I(t)')
    plt.grid(True)
    plt.legend()

plt.show()
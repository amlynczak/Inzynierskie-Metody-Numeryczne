import numpy as np
import matplotlib.pyplot as plt

from  parameters import beta, N, gamma, tmax, dt, u0, TOL, max_iterations, f

# Metoda Newtona do rozwiązywania równań predyktora
def newton_predictor(U1, U2, un, dt, alpha, beta, max_iterations, TOL):
    for _ in range(max_iterations):
        F1 = U1 - un - dt * (0.5 * alpha * U1 - 0.5 * beta * U2) #wzor (18)
        F2 = U2 - un - dt * (0.5 * alpha * U2 - 0.5 * beta * U2) #wzor (19)
        
        m1_1 = 1 - 0.5 * dt * alpha
        m1_2 = -0.5 * dt * beta
        m2_1 = -0.5 * dt * alpha
        m2_2 = 1 - 0.5 * dt * beta 
        determinant = m1_1 * m2_2 - m1_2 * m2_1 #cala sekwencja wylicza wedlug wzorów (22) - (26)
        
        delta_U1 = (F2 * m1_2 - F1 * m2_2) / determinant
        delta_U2 = (F1 * m2_1 - F2 * m1_1) / determinant #wzory (27) oraz (28)
        
        U1 = U1 + delta_U1
        U2 = U2 + delta_U2
        
        if abs(delta_U1) < TOL and abs(delta_U2) < TOL: #warunek stopu 
            break
    
    return U1, U2

# Dwuetapowa metoda RK2
def rk2_method(u0, tmax, dt, TOL, max_iterations):
    t = np.arange(0, tmax + dt, dt) #rowne rozlozenie t w zakresie
    u = [u0] #wektor dla u(t)
    
    for i in range(len(t) - 1):
        un = u[-1]
        alpha = beta * N - gamma
        
        # Równania predyktora
        U1 = un + dt * (0.5 * f(t[i], un) + 0.5 * f(t[i] + 0.5 * dt, un))
        U2 = un + dt * (0.5 * f(t[i], un) + 0.5 * f(t[i] + 0.5 * dt, U1))
        
        U1, U2 = newton_predictor(U1, U2, un, dt, alpha, beta, max_iterations, TOL)
        
        unew = un + dt * (0.5 * f(t[i], un) + 0.5 * f(t[i] + dt, U2)) # Obliczenie u_n+1, wzor (17)
        u.append(unew)
    
    return t, u

t_rk2, u_rk2 = rk2_method(u0, tmax, dt, TOL, max_iterations) #uzycie metody
z_rk2 = [N - u_t for u_t in u_rk2] #obliczenie z(t) na podstawie otrzymanego u(t)

plt.figure(figsize=(10, 6))
plt.plot(t_rk2, u_rk2, label="u(t)")
plt.plot(t_rk2, z_rk2, label="z(t)")
plt.xlabel("t")
plt.ylabel("u(t), z(t)")
plt.legend()
plt.title("Niejawna RK2")
plt.grid()
plt.show()
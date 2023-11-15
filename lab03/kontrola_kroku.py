import matplotlib.pyplot as plt
from plots import plot_results
from rk2 import rk2
from trapez import metoda_trapezow
from fundamentals import stala_bledu
from fundamentals import x0, v0, dt0
from fundamentals import tmax, alpha, S, p, TOL_values

#funckja zmiany kroku
def zmiana_kroku(S, TOL, max, p, dt):
    tmp = (S*TOL)/max
    ret = tmp**(1/(p+1)) * dt
    return ret

#glowna funkcja kontroli kroku - calosc wedlug pseudokodu dostarczonego w instrukcji
def kontrola_kroku(x0, v0, dt0, tmax, TOL, S, p, fun_num):
    t = 0
    dt = dt0
    xn = x0
    vn = v0

    t_values = []
    dt_values = []
    xn_values = []
    vn_values = []

    while t < tmax:
        x_2n1, v_2n1 = fun_num(xn, vn, t, alpha, dt, TOL)
        x_2n2, v_2n2 = fun_num(x_2n1, v_2n1, t, alpha, dt, TOL)

        x_1n2, v_1n2 = fun_num(xn, vn, t, alpha, 2*dt, TOL)

        Ex = stala_bledu(x_1n2, x_2n2, p)
        Ev = stala_bledu(v_1n2, v_2n2, p)

        if max(abs(Ex), abs(Ev)) < TOL:
            t = t + 2*dt
            xn = x_2n2
            vn = v_2n2

            t_values.append(t)
            dt_values.append(dt)
            xn_values.append(xn)
            vn_values.append(vn)
        
        dt = zmiana_kroku(S, TOL, max(abs(Ex), abs(Ev)), p, dt)

    return t_values, dt_values, xn_values, vn_values

# w jednej funkcji wywolanie kontroli korku oraz wyswietlenia wykresow
def rozwiazanie_z_konrola(x0, v0, dt0, tmax, alpha, TOL_values, fun_num):
    TOL = TOL_values[0]
    t_values_1, dt_values_1, x_values_1, v_values_1 = kontrola_kroku(x0, v0, dt0, tmax, TOL, S, p, fun_num)
    TOL = TOL_values[1]
    t_values_2, dt_values_2, x_values_2, v_values_2 = kontrola_kroku(x0, v0, dt0, tmax, TOL, S, p, fun_num)

    if fun_num == metoda_trapezow:
        title = 'Metoda Trapezow'
    if fun_num == rk2:
        title = 'Metoda RK2'
    
    plot_results(t_values_1, x_values_1, v_values_1, dt_values_1, t_values_2, x_values_2, v_values_2, dt_values_2, title)


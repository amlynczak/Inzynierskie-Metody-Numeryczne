from fundamentals import f, g
import numpy as np

#implementacja schematu numerycznego Metody Trapezow
def metoda_trapezow(xn, vn, tn, alpha, dt, TOL):
    # definicja dwoch funkcji nieliniowych F i G - wzory (11) i (12)
    F = lambda x_n1, v_n1: x_n1 - xn - (dt/2) * (f(tn, xn, vn) + f(tn+1, x_n1, v_n1))
    G = lambda x_n1, v_n1: v_n1 - vn - (dt/2) * (g(tn, xn, vn, alpha) + g(tn+1, x_n1, v_n1, alpha))

    # inicjalizaca wartosci, ktore bedzemy zwracac (przyjmujemy je startowo jako rowne xn i vn)
    x_0n1 = xn
    v_0n1 = vn

    # iteracje
    while True:
        # Elementy macierzy ukladu A - wzory (16) - (19)
        a11 = 1
        a12 = -dt/2
        a21 = -dt/2 * (-2 * alpha * x_0n1 * v_0n1 - 1)
        a22 = 1 - dt/2 * alpha * (1 - x_0n1**2)

        # liczymy wyznacznik macierzy do wzorow (21) oraz (22) - mianownik w jednym i drugim
        det = a11 * a22 - a12 * a21

        # obliczanie dt oraz dv
        delta_x = ((-F(x_0n1, v_0n1)) * a22 - (-G(x_0n1, v_0n1)) * a12) / det #wzor (21)
        delta_v = (a11 * (-G(x_0n1, v_0n1)) - a21 * (-F(x_0n1, v_0n1))) / det #wzor (22)

        # Update 
        x_0n1 += delta_x
        v_0n1 += delta_v

        # sprawdzamy warunek konca petli
        if abs(delta_x) < TOL and abs(delta_v) < TOL:
            break

    return x_0n1, v_0n1
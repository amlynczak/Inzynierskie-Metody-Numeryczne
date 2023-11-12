# Parametry startowe
x0 = 0.01
v0 = 0
dt0 = 1
tmax = 40
alpha = 5
S = 0.75
p = 2

# Parametry dla kontrolowanej kroku - dwie wartosci
TOL_values = [1e-2, 1e-5]

#funkcja f - wzor (2)
def f(t, x, v):
    return v

#funckja g - wzor (3)
def g(t, x, v, alpha):
    return alpha * (1 - x**2) * v - x

#potrzebne do implementacji i obliczenia Ex, Ev
def stala_bledu(x1, x2, p):
    err = (x2 - x1)/(2**p - 1)
    return err
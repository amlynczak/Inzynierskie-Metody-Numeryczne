# Parametry
beta = 0.001
N = 500
gamma = 0.1
tmax = 100
dt = 0.1
u0 = 1
TOL = 1e-6
max_iterations = 20

# Funkcja f(t, u), wzór (2)
def f(t, u):
    alpha = beta * N - gamma
    return (alpha * u - beta * u**2)
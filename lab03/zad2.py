from kontrola_kroku import rozwiazanie_z_konrola
from fundamentals import x0, v0, dt0, tmax, alpha, TOL_values
from trapez import metoda_trapezow

rozwiazanie_z_konrola(x0, v0, dt0, tmax, alpha, TOL_values, metoda_trapezow)
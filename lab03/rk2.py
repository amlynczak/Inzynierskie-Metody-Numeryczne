from fundamentals import f, g

#implementacja schematu numerycznego metody RK2 - instukcja wzory (25) do (31)
def rk2(xn, vn, tn, alpha, dt, TOL):
    k1_x = f(tn, xn, vn)
    k1_v = g(tn, xn, vn, alpha)

    k2_x = f(tn + dt, xn + dt*k1_x, vn + dt * k1_v)
    k2_v = g(tn + dt, xn + dt*k1_x, vn + dt * k1_v, alpha)

    xn1 = xn + (dt/2)*(k1_x + k2_x)
    vn1 = vn + (dt/2)*(k1_v + k2_v)

    return xn1, vn1
import numpy as np
from rungekutta4 import rk4
from euler import euler


def dydx(x,y):
    return x**2 - 3*y + 3*x*y**2

for i in range(1,10):
    p = 10**(-i)
    r = rk4(dydx, 0, 2, 0, p)
    e = euler(dydx, 0, 2, 0, p)
    err = abs(r - e)
    print(f"p = {p}, r = {r}, e = {e}, err = {err}")
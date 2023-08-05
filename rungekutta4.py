import numpy as np
import matplotlib.pyplot as plt

def rk4(func, x0, xf,y0, h):
    x = x0
    y = y0
    i = 0
    while x < xf:
        k1 = func(x,y)
        k2 = func(x + h/2, y + h*k1/2)
        k3 = func(x + h/2, y + h*k2/2)
        k4 = func(x + h, y + h*k3)
        y = (y + h*(k1 + 2*k2 + 2*k3 + k4)/6)
        x = x + h
        i = i + 1
        #$print(f"x: {x}, k1: {k1}, k2: {k2}, k3: {k3}, k4: {k4}")
        #print(f"i = {i}, x = {x}, y = {y}")
    return y

# def dydx(x,y):
#     return 2*x + 3*y

# for i in range (1,5):
#     h = 10**(-i)
#     y = rk4(dydx, 0, 2, 0, h)
#     print(y)


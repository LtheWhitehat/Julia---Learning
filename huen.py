import numpy as np
from scipy.integrate import solve_ivp
from euler import euler

def euler_system(f, t0, y0, h, N):
    t_values = [t0]
    y_values = [y0]

    for i in range(N):
        t_i = t_values[-1]
        y_i = y_values[-1]
        y_i_plus_1 = [y_i[k] + h * f(t_i, y_i)[k] for k in range(2)]

        t_values.append(t_i + h)
        y_values.append(y_i_plus_1)

    return t_values, y_values


def heun(dydx,x0,y0,xf,h):
    n = int((xf - x0) / h)
    print(n)
    y = [y0]
    x = [x0]
    for i in range(1,n) :
        k1 = dydx(x[-1],y[-1])
        ye = y[-1] + h*k1
        k2 = dydx(x[-1], ye)
        y.append(y[-1] + (h/2)*(k1 + k2))
        x.append(x[-1] + h)
        #print(x0,y0)
        
    return x,y

def rk2(dydx,x0,y0,xf,h):
    n = int((xf - x0) / h)
    x = [x0]
    y = [y0]
    for i in range(1,n):
        k1 = dydx(x[-1],y[-1])*h
        k2 = h*dydx(x[-1] + h, y[-1] + k1*h)
        y.append(y[-1] + (k1 + 2*k2)/6)
        x.append(x[-1] + h)
        print(f"x: {x[-1]}, y: {y[-1]}")
        
    return x,y


def f(y, v, x):
    return y+v+3*x
def g(y, v, x):
    return 2*y-v-x

def euler_system(y0, v0, x0, h, num_steps):
    y = y0
    v = v0
    x = x0
    for _ in range(num_steps):
        y_new = y + h * f(y, v, x)
        v_new = v + h * g(y, v, x)
        y = y_new
        v = v_new
        x += h
        print(f"x: {x}, y: {y}, v: {v}")
    return y, v

# Initial conditions and parameters
y0 = 0.
v0 = -1.
x0 = 0.
h = 0.2  # Time step size
num_steps = 2  # Number of steps

# Solve the second-order ODE using Euler's method
final_y, final_v = euler_system(y0, v0, x0, h, num_steps)

print("Final y:", final_y)

def dydx(x,y):
    return -3*x**2*y

x0 = 0
y0 = 2
h = 0.5

euler(dydx, x0, 1, y0, h)



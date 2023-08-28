from metodos_corretor_preditor import *
from rungekutta4 import rk4
from euler import euler
from scipy.integrate import solve_ivp
import numpy as np

t0 = 0
y0 = 0  # Initial values for y and z
u0 = 1
h = 0.1
N = 10

def odes_system(t, y):
    u, v = y[0], y[1]
    print(u, v)
    du_dt = v
    dv_dt = -u * v - u
    return [du_dt, dv_dt]



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


def odes_system1(t, y):
    y_prime = y[0] + y[1] + 3 * t
    z_prime = 2 * y[0] - y[1] - t
    return [y_prime, z_prime] # Return individual derivative values, not a list

t0 = 0  # Initial time
y0 = [0, -1]  # Initial values for y, y', z, z'
h = 0.2  # Step size
N = 6  # Number of steps


import numpy as np
from scipy.integrate import solve_ivp

def system(t, y):
    """
    Este código implementa o sistema de EDOs:
    y' = y + z + 3*x
    z' = 2y - z - x
    """
    dy = y + z + 3*x
    dz = 2*y - z - x
    return dy, dz

import numpy as np

# Define the system of ODEs as functions
def f1(x, y, z):
    return  y + z + 3*t

def f2(x, y, z):
    return 2*y - z - x

# Implement the RK4 method for a single time step
def rk4_step(t, dt, y1, y2):
    k1_y1 = dt * f1(t, y1, y2)
    k1_y2 = dt * f2(t, y1, y2)
    k2_y1 = dt * f1(t + 0.5*dt, y1 + 0.5*k1_y1, y2 + 0.5*k1_y2)
    k2_y2 = dt * f2(t + 0.5*dt, y1 + 0.5*k1_y1, y2 + 0.5*k1_y2)
    k3_y1 = dt * f1(t + 0.5*dt, y1 + 0.5*k2_y1, y2 + 0.5*k2_y2)
    k3_y2 = dt * f2(t + 0.5*dt, y1 + 0.5*k2_y1, y2 + 0.5*k2_y2)
    k4_y1 = dt * f1(t + dt, y1 + k3_y1, y2 + k3_y2)
    k4_y2 = dt * f2(t + dt, y1 + k3_y1, y2 + k3_y2)

    y1_new = y1 + (1/6) * (k1_y1 + 2*k2_y1 + 2*k3_y1 + k4_y1)
    y2_new = y2 + (1/6) * (k1_y2 + 2*k2_y2 + 2*k3_y2 + k4_y2)

    print(f"t: {t}, y1: {y1_new}, y2: {y2_new}")
    return y1_new, y2_new

# Set up initial conditions
t0 = 0.0
y1_0 = 0.0
y2_0 = -1.0

# Choose the integration time step and range
dt = 0.2
t_end = 2

# Create lists to store the time and solution values
time_points = [t0]
y1_values = [y1_0]
y2_values = [y2_0]

# Perform the integration using RK4
t = t0
y1 = y1_0
y2 = y2_0
while t < t_end:
    y1, y2 = rk4_step(t, dt, y1, y2)
    t += dt

    # Append the results to the lists
    time_points.append(t)
    y1_values.append(y1)
    y2_values.append(y2)



def system_of_odes(t, y):
    x, y = y
    dxdt = 2 * x + 4 * y
    dydt = -x + 6 * y
    return [dxdt, dydt]


def dydx(x,y):
    return -3*x**2*y
    # return -1 + y/x
    #return x**3 + x**2 + x + 1
    #return y*np.log(y)/x
x0 = 0
y0 = 2
h = 0.5

#print(euler(dydx, x0, 1, y0, h))



# sol = solve_ivp(dydx, (t0, t0 + h), y0, method='RK23',dense_output=  True)
# print(sol.sol(0.1))

# t,y,z = (adams_bashforth_multistep_rk4(dydx, x0, y0, h, 2))
# for t,y in zip(t,y):
#     print(f"t: {t}, y: {y}")
# x,y_pred,y = adams_bashforth_multistep_rk4(dydx, x0, y0, h,6)

# for x,y_pred,y in zip(x,y_pred,y):
#     print(f"x: {x}, y_pred: {y_pred}, y: {y}")

#print(euler(dydx, x0, 2.3, y0, h))

# # Now you have the solutions in time_points, y1_values, and y2_values
# for x,y,z in zip(time_points, y1_values, y2_values):
#     print(f"t = {x:.2f}, y1 = {y:.6f}, y2 = {z:.6f}")



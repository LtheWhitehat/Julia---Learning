from metodos_corretor_preditor import *
from rungekutta4 import rk4
from euler import euler
from scipy.integrate import solve_ivp


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

# t_values = [t0]
# y_values = [y0, u0]

# for i in range(N):
#     t_i = t_values[-1]
#     y_i = y_values[-1]
#     y_i_plus_1 = rk4_system(odes_system, t_i, y_i, h)
#     t_values.append(t_i + h)
#     y_values.append(y_i_plus_1)
    
    
# for t, y in zip(t_values, y_values):
#     print(f"t: {t}, y: {y}")
    
    
#print(euler(f, t0, 2, y0, h))

# x_val,y_val = adams_moulton_rk4(f,t0, y0, h, N)
# x_bash, y_bash = adams_bashforth_multistep(f,t0,y0,h,N)
# for xm, ym, xb, yb in zip(x_val, y_val, x_bash, y_bash):
#     print(f"Bash x: {xb}, y: {yb}\n")
#     #print(f"Moulton x: {xm} -- y: {ym}\nBash x: {xb}, y: {yb}\n\n")


# def rk4_system(f, t, y, h):
#     k1 = h * f(t, y)
#     k2 = h * f(t + h / 2, y + k1 / 2)
#     k3 = h * f(t + h / 2, y + k2 / 2)
#     k4 = h * f(t + h, y + k3)
#     return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6

# def f(t, y):
#     dydt = y[1]
#     dzdt = -y[0]
#     return [dydt, dzdt]



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

sol = solve_ivp(odes_system1, [t0, t0 + h * N], y0, t_eval=np.linspace(t0, t0 + h * N, N + 1), dense_output=True)
#sol.sol(0.2)
print(sol.sol(0.2))

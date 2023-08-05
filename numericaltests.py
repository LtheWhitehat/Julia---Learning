from metodos_corretor_preditor import *
from rungekutta4 import rk4
from euler import euler


def f(x, y):
    return np.sin(x+y) - np.exp(x)

t0 = 0  # Initial time
y0 = 4  # Initial value
h = 0.5 # Step size
N = 5  # Number of steps

print(euler(f, t0, 2, y0, h))

# x_val,y_val = adams_moulton_rk4(f,t0, y0, h, N)
# x_bash, y_bash = adams_bashforth_multistep(f,t0,y0,h,N)
# for xm, ym, xb, yb in zip(x_val, y_val, x_bash, y_bash):
#     print(f"Bash x: {xb}, y: {yb}\n")
#     #print(f"Moulton x: {xm} -- y: {ym}\nBash x: {xb}, y: {yb}\n\n")


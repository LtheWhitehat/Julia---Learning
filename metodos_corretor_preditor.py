import numpy as np
from euler import euler
#from rungekutta4 import rk4

def adams_bashforth_multistep(f, t0, y0, h, N):
    t_values = [t0]
    y_values = [y0]

    # Use an initial method (e.g., Euler's method) to get the second point
    t_values.append(t0 + h)
    y_values.append(y_values[-1] + h*f(t_values[-1], y_values[-1]))
    y_pred = y_values
    y_corr = y_values
    for i in range(2, N + 1):
        t_i = t_values[-1]
        y_i = y_values[-1]

        # Adams-Bashforth predictor
        predictor_sum = 0
        for j in range(1, i + 1):
            predictor_sum += 3 * f(t_values[-j], y_values[-j])
        y_predictor = y_i + (h / 2) * predictor_sum
        y_pred.append(y_predictor)
        # Adams-Moulton corrector (implicit formula)
        y_i_plus_1 = y_i + (h / 2) * (f(t_i, y_i) + f(t_i + h, y_predictor))
        y_corr.append(y_i_plus_1)
        t_values.append(t_i + h)
        y_values.append(y_i_plus_1)

    return t_values, y_values

def rk4(f, t, y, h):
    k1 = h * f(t, y)
    k2 = h * f(t + h / 2, y + k1 / 2)
    k3 = h * f(t + h / 2, y + k2 / 2)
    k4 = h * f(t + h, y + k3)
    return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6


def adams_bashforth_multistep_rk4(f, t0, y0, h, N):
    t_values = [t0]
    y_values = [y0]

    # Use RK4 to get the second point
    t_values.append(t0 + h)
    y_values.append(rk4(f, t0, y0, h))
    y_pred = y_values
    for i in range(2, N + 1):
        t_i = t_values[-1]
        y_i = y_values[-1]

        # Adams-Bashforth predictor
        predictor_sum = 0
        for j in range(1, min(i, 4) + 1):  # Use up to the 4th order predictor (Adams-Bashforth)
            predictor_sum += (3 if j == 1 else 2) * f(t_values[-j], y_values[-j])

        y_predictor = y_i + (h / 2) * predictor_sum
        y_pred.append(y_predictor)
        # Adams-Moulton corrector (implicit formula)
        y_i_plus_1 = y_i + (h / 2) * (f(t_i + h, y_i) + f(t_i + h, y_predictor))

        t_values.append(t_i + h)
        y_values.append(y_i_plus_1)

    return t_values, y_values, y_pred


def adams_moulton_rk4(f, t0, y0, h, N):
    t_values = [t0]
    y_values = [y0]

    # Use RK4 to get the second point
    t_values.append(t0 + h)
    y_values.append(rk4(f, t0, y0, h))

    for i in range(2, N + 1):
        t_i = t_values[-1]
        y_i = y_values[-1]

        # Adams-Moulton corrector (implicit formula)
        predictor_sum = 0
        for j in range(1, i + 1):
            predictor_sum += f(t_values[-j], y_values[-j])
        y_i_plus_1 = y_i + (h / 2) * (f(t_i + h, y_i) + predictor_sum)

        t_values.append(t_i + h)
        y_values.append(y_i_plus_1)

    return t_values, y_values

def adams_bashforth_2(f, t0, y0, h, N):
    t_values = [t0]
    y_values = [y0]

    for i in range(1,N):
        t_i = t_values[-1]
        y_i = y_values[-1]
        y_bash = []
        y_corr = []
        # Adams-Bashforth predictor
        f_predictor = f(t_i, y_i)
        y_predictor = y_i + h * f_predictor
        y_bash.append(y_predictor)
        # Adams-Moulton corrector (implicit formula)
        f_corrector = f(t_i + h, y_predictor)
        y_corr.append(f_corrector)
        y_i_plus_1 = y_i + h * (f_predictor + f_corrector) / 2

        t_values.append(t_i + h)
        y_values.append(y_i_plus_1)

    return t_values, y_values, y_bash, y_corr

# Example usage:
# Define the ODE: dy/dt = f(t, y)

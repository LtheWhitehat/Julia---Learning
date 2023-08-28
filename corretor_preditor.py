import numpy as np

def f(t, y):
    """
    Define a função que representa a equação diferencial.

    Parameters:
        t : float
            O tempo atual.
        y : float
            O valor atual da variável y.

    Returns:
        dydt : float
            O valor da derivada de y no tempo t.
    """
    return -t*y

def adams_bashforth_4_step(t, dt, y):
    """
    Realiza um passo do método de Adams-Bashforth de passo 4.

    Parameters:
        t : float
            O tempo atual.
        dt : float
            O tamanho do passo de tempo.
        y : float
            O valor atual da variável y.

    Returns:
        y_new : float
            O valor de y no tempo t + dt.
    """
    f0 = f(t, y)
    f1 = f(t - dt, y - dt * f0)
    f2 = f(t - 2*dt, y - 2*dt * f1)
    f3 = f(t - 3*dt, y - 3*dt * f2)
    y_new = y + dt/24 * (55*f0 - 59*f1 + 37*f2 - 9*f3)
    return y_new

def adams_moulton_3_step(t, dt, y):
    """
    Realiza um passo do método de Adams-Moulton de passo 3.

    Parameters:
        t : float
            O tempo atual.
        dt : float
            O tamanho do passo de tempo.
        y : float
            O valor atual da variável y.

    Returns:
        y_new : float
            O valor de y no tempo t + dt.
    """
    f0 = f(t, y)
    f1 = f(t - dt, y - dt * f0)
    f2 = f(t - 2*dt, y - 2*dt * f1)
    y_new = y + dt/24 * (9*f(t + dt, adams_bashforth_4_step(t, dt, y)) + 19*f0 - 5*f1 + f2)
    return y_new

# Set up initial conditions
t0 = 0.0
y0 = 1.0

# Choose the integration time step and range
dt = 0.5
t_end = 2.5

# Create lists to store the time and solution values
time_points = [t0]
solution = [y0]

# Perform the integration using Adams-Bashforth of order 4 for the first step
t = t0
y = y0
while t < t0 + 3*dt:  # Using Adams-Bashforth for the first 3 steps
    y = adams_bashforth_4_step(t, dt, y)
    t += dt

    # Append the results to the lists
    time_points.append(t)
    solution.append(y)

# Continue the integration using Adams-Moulton of order 3
while t < t_end:
    y = adams_moulton_3_step(t, dt, y)
    t += dt

    # Append the results to the lists
    time_points.append(t)
    solution.append(y)
    print(f"t = {t:.2f}, y = {y:.6f}")

# Now you have the values of y at different time points

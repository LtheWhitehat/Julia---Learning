import numpy as np

def system_of_odes(t, y):
    """
    Define the system of ODEs.

    Parameters:
        t : float
            The current time.
        y : array-like
            An array-like object containing the state variables x and y.

    Returns:
        dydt : array-like
            An array-like object containing the derivatives of x and y at time t.
    """
    x, y = y
    dxdt = 2 * x + 4 * y
    dydt = -x + 6 * y
    return [dxdt, dydt]

def euler_step(t, dt, y):
    dydt = np.array(system_of_odes(t, y))
    y_new = y + dt * dydt
    return y_new

# Set up initial conditions
t0 = 0.0
y0 = [-1, 6]

# Choose the integration time step and range
dt = 0.1
t_end = 1.0

# Create lists to store the time and solution values
time_points = [t0]
solution = [y0]

# Perform the integration using Euler's method
t = t0
y = y0
while t < t_end:
    y = euler_step(t, dt, y)
    t += dt

    # Append the results to the lists
    time_points.append(t)
    solution.append(y)
    print(t,y)



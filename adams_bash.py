def f(t, y):
    """
    Define the function representing the ODE dy/dt = f(t, y).

    Parameters:
        t : float
            The current time.
        y : float
            The current value of the solution.

    Returns:
        dydt : float
            The value of the derivative dy/dt at time t.
    """
    return t+y-1

def adams_bashforth(t, dt, y):
    """
    Perform one step of the Adams-Bashforth method.

    Parameters:
        t : float
            The current time.
        dt : float
            The time step size.
        y : float
            The current value of the solution.

    Returns:
        y_new : float
            The value of the solution at time t + dt.
    """
    print(f"t = {t}, y = {y}")
    f_n = f(t, y)
    f_n_minus_1 = f(t - dt, y - dt * f_n)
    y_new = y + dt * (0.5 * f_n + 0.5 * f_n_minus_1)
    return y_new

# Set up initial conditions
t0 = 0
y0 = 1

# Choose the integration time step and range
dt = 0.2
t_end = 1.2

# Create lists to store the time and solution values
time_points = [t0]
solution = [y0]

# Perform the integration using Adams-Bashforth method
t = t0
y = y0
while t < t_end:
    y = adams_bashforth(t, dt, y)
    t += dt

    # Append the results to the lists
    time_points.append(t)
    solution.append(y)

# Now you have the values of the solution at different time points
for val in solution:
    print(val)
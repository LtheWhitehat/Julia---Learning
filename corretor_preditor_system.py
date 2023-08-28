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

def adams_bashforth_4_step(t, dt, y_history):
    """
    Perform one step of the Adams-Bashforth method of order 4.

    Parameters:
        t : float
            The current time.
        dt : float
            The time step size.
        y_history : list
            A list containing the history of state variables at previous time points.

    Returns:
        y_new : array-like
            An array-like object containing the updated state variables x and y at time t + dt.
    """
    f_t = np.array(system_of_odes(t, y_history[-1]))
    f_t_minus_1 = np.array(system_of_odes(t - dt, y_history[-2]))
    f_t_minus_2 = np.array(system_of_odes(t - 2*dt, y_history[-3]))
    f_t_minus_3 = np.array(system_of_odes(t - 3*dt, y_history[-4]))

    y_new = y_history[-1] + (dt/24) * (55*f_t - 59*f_t_minus_1 + 37*f_t_minus_2 - 9*f_t_minus_3)
    return y_new

def adams_moulton_3_step(t, dt, y_history):
    """
    Perform one step of the Adams-Moulton method of order 3.

    Parameters:
        t : float
            The current time.
        dt : float
            The time step size.
        y_history : list
            A list containing the history of state variables at previous time points.

    Returns:
        y_new : array-like
            An array-like object containing the updated state variables x and y at time t + dt.
    """
    f_t_plus_1_predicted = np.array(system_of_odes(t + dt, adams_bashforth_4_step(t, dt, y_history[-4:])))
    f_t = np.array(system_of_odes(t, y_history[-1]))
    f_t_minus_1 = np.array(system_of_odes(t - dt, y_history[-2]))

    y_new = y_history[-1] + (dt/24) * (9*f_t_plus_1_predicted + 19*f_t - 5*f_t_minus_1)
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

# Perform the integration using Adams-Bashforth of step 4 and Adams-Moulton of step 3
t = t0
y = y0
while t < t_end:
    if len(solution) < 4:
        y = adams_bashforth_4_step(t, dt, solution)
    else:
        y = adams_moulton_3_step(t, dt, solution)
    t += dt

    # Append the results to the lists
    time_points.append(t)
    solution.append(y)

# Extract x and y values from the solution list
x_values, y_values = zip(*solution)

# Now you have the values of x and y at different time points

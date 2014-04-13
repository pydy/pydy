#!/usr/bin/env python

# This is an example of integrating the equations of motion for a double
# pendulum which were generated with sympy.physics.mechanics. We make use of
# SciPy/NumPy for the integration routines and Matplotlib for plotting.
#
# Steps taken:
# 1. Turned on mechanics_printing() in sympy.physics.mechanics for proper
# output for copying the equations to this file.
# 2. Import zeros, sin, cos, linspace from NumPy, odeint from SciPy and pyplot
# from Matplotlib.
# 3. Write a function definition that returns the right hand side of the
# first order form of the equations of motion. rhd(y, t, *parameters)
# 4. Called odeint with rhs, y0, t and the parameters.
# 5. Plotted the results.

from numpy import sin, cos, linspace, zeros
import matplotlib.pyplot as plt
from scipy.integrate import odeint

## Integration ##

def rhs(y, t, l, m, g):
    """Returns the derivatives of the states at the given time for the given
    set of parameters.

    Parameters
    ----------
    y : array_like, shape(n,)
        An array of the current states.
    t : float
        The current time.
    l : float
        Pendulum length.
    m : float
        Pendulum mass.
    g : float
        Acceleration due to gravity.

    Returns
    -------
    dydt : array_like, shape(n,)
        An array of the current derivatives of the states.

    Notes
    -----
    The units and order of the states, time and parameters should be
    consistent.

    """
    # Unpack the states so you can use the variable names in the
    # sympy.physics.mechanics equations
    q1 = y[0]
    q2 = y[1]
    u1 = y[2]
    u2 = y[3]
    # or you can make use of python's tuple unpacking for a one liner
    # q1, q2, u1, u2 = y

    # Initialize a vector for the derivatives.
    dydt = zeros((len(y)))

    # Compute the derivatives, these are pasted in from the
    # sympy.physics.mechanics results.
    dydt[0] = u1
    dydt[1] = u2
    dydt[2] = (-g*sin(q1)*sin(q2)**2 + 2*g*sin(q1) -
        g*sin(q2)*cos(q1)*cos(q2) + 2*l*u1**2*sin(q1)*cos(q1)*cos(q2)**2 -
        l*u1**2*sin(q1)*cos(q1) - 2*l*u1**2*sin(q2)*cos(q1)**2*cos(q2) +
        l*u1**2*sin(q2)*cos(q2) + l*u2**2*sin(q1)*cos(q2) -
        l*u2**2*sin(q2)*cos(q1))/(l*(sin(q1)**2*sin(q2)**2 +
        2*sin(q1)*sin(q2)*cos(q1)*cos(q2) + cos(q1)**2*cos(q2)**2 - 2))
    dydt[3] = (-sin(q1)*sin(q2)/2 - cos(q1)*cos(q2)/2)*(2*g*l*m*sin(q1) -
        l**2*m*(-sin(q1)*cos(q2) +
        sin(q2)*cos(q1))*u2**2)/(l**2*m*(sin(q1)*sin(q2)/2 +
        cos(q1)*cos(q2)/2)*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) -
        l**2*m) + (g*l*m*sin(q2) - l**2*m*(sin(q1)*cos(q2) -
        sin(q2)*cos(q1))*u1**2)/(l**2*m*(sin(q1)*sin(q2)/2 +
        cos(q1)*cos(q2)/2)*(sin(q1)*sin(q2) + cos(q1)*cos(q2))
        - l**2*m)

    # Return the derivatives.
    return dydt

# Specify the length, mass and acceleration due to gravity.
parameters = (1, 1, 9.8)
# Specify initial conditions for the states.
y0 = [.1, .2, 0, 0]
# Create a time vector.
t = linspace(0, 5)
# Integrate the equations of motion.
y = odeint(rhs, y0, t, parameters)


## Plotting ##

# Create an empty figure.
fig = plt.figure()

# Add a single axes to the figure.
ax = fig.add_subplot(1, 1, 1)

# Plot the states versus time.
ax.plot(t, y)

# Add a title, axes labels and a legend.
ax.set_title('Double Pendulum Example')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Angle, Angluar rate (rad, rad/s)')
ax.legend(['q1', 'q2', 'u1', 'u2'])

# Display the figure.
plt.show()

# turned on mechanics_printing()
# Steps: wrote function definition; f(y, t, arguements)
# Wrote return statement, copy/pasted from ipython console and formatted it
# unpacked y
# had to import array, sin, cos, linspace from numpy, odeint from scipy
# added args, y0, t vector, odeint call
# added plotting code

from numpy import array, sin, cos, linspace
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def rhs(y, t, l, m, g):
    q1, q2, u1, u2 = y
    return array([u1, u2, (-g*sin(q1)*sin(q2)**2 + 2*g*sin(q1) -
        g*sin(q2)*cos(q1)*cos(q2) + 2*l*u1**2*sin(q1)*cos(q1)*cos(q2)**2 -
        l*u1**2*sin(q1)*cos(q1) - 2*l*u1**2*sin(q2)*cos(q1)**2*cos(q2) +
        l*u1**2*sin(q2)*cos(q2) + l*u2**2*sin(q1)*cos(q2) -
        l*u2**2*sin(q2)*cos(q1))/(l*(sin(q1)**2*sin(q2)**2 +
        2*sin(q1)*sin(q2)*cos(q1)*cos(q2) + cos(q1)**2*cos(q2)**2 - 2)),
        (-sin(q1)*sin(q2)/2 - cos(q1)*cos(q2)/2)*(2*g*l*m*sin(q1) -
        l**2*m*(-sin(q1)*cos(q2) +
        sin(q2)*cos(q1))*u2**2)/(l**2*m*(sin(q1)*sin(q2)/2 +
        cos(q1)*cos(q2)/2)*(sin(q1)*sin(q2) + cos(q1)*cos(q2)) -
        l**2*m) + (g*l*m*sin(q2) - l**2*m*(sin(q1)*cos(q2) -
        sin(q2)*cos(q1))*u1**2)/(l**2*m*(sin(q1)*sin(q2)/2 +
        cos(q1)*cos(q2)/2)*(sin(q1)*sin(q2) + cos(q1)*cos(q2))
        - l**2*m)])

# args are l, m, g
args = (1, 1, 9.8)
y0 = [.1, .2, 0, 0]
t = linspace(0, 5)
y = odeint(rhs, y0, t, args)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t, y)
ax.set_title('Double Pendulum Example')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Angle, Angluar rate (rad, rad/s)')
ax.legend(['q1', 'q2', 'u1', 'u2'])
plt.show()

"""
# How to integrate equations of motion, quick and dirty way
# get EoM into form <qdots, udots = expressions> first though
# Also, make sure there are no qdots in rhs of udots
# (meaning udot = f(q, u, t), not f(q, qdot, u, t)
# use Kane.kindiffdict to get dictionary, and use subs on your udot vector to
# get rid of qdots in bad places. See examples

from numpy import array, linspace, sin, cos # etc.
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def rhs(y, t, arg1, arg2, etc.):
    # unpack variables; remember the order here, it is important
    var1, var2, etc. = y
    # write return statement
    # copy/paste qdots, udots from console/terminal
    # make sure mechanics printing was turned on first though
    return array([var1dot, var2dot, var3dot, etc.])

# make tuple of arguement values (as floats or something)
args = (arg1, arg2, etc.)
# give initial conditions
y0 = []
# choose a timespan
t = linspace(0, 1)
# call odeint
y = odeint(rhs, y0, t, args)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t, y)
# give title to plot
ax.set_title('Double Pendulum Example')
# give x axis label
ax.set_xlabel('Time (s)')
# give y axis label
ax.set_ylabel('Angle, Angluar rate (rad, rad/s)')
# set legend values
ax.legend(['var1', 'var2', 'var3', etc.])
# show plot
plt.show()
"""

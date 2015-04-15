# How to integrate equations of motion, quick and dirty way
# Note: this template will not run as is
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
ax.set_title('My Problem')
# give x axis label
ax.set_xlabel('Time (s)')
# give y axis label
ax.set_ylabel('Qs and Us')
# set legend values
ax.legend(['var1', 'var2', 'var3', etc.])
# show plot
plt.show()

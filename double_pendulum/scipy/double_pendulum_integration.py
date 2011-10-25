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



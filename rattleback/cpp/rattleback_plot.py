import numpy as np
import matplotlib.pyplot as plt
#import os

# Define a data type for the simdata structure
simdata = np.dtype([('t', np.float64),
                    ('q0', np.float64),
                    ('q1', np.float64),
                    ('q2', np.float64),
                    ('q3', np.float64),
                    ('q4', np.float64),
                    ('u0', np.float64),
                    ('u1', np.float64),
                    ('u2', np.float64),
                    ('qd0', np.float64),
                    ('qd1', np.float64),
                    ('qd2', np.float64),
                    ('qd3', np.float64),
                    ('qd4', np.float64),
                    ('ud0', np.float64),
                    ('ud1', np.float64),
                    ('ud2', np.float64),
                    ('Rx', np.float64),
                    ('Ry', np.float64),
                    ('Rz', np.float64),
                    ('ke', np.float64),
                    ('pe', np.float64),
                    ('te', np.float64),
                    ('delta', np.float64)])

#os.system('make simulation.dat')
data = np.fromfile('simulation.dat', dtype=simdata) # read the data

plt.figure()
plt.subplot(211)
plt.plot(data['t'], data['delta']*180./np.pi)
plt.ylabel(r'$\delta$')
plt.subplot(212)
plt.plot(data['t'], data['q0']*180./np.pi)
plt.ylabel(r'$\gamma$')
plt.xlabel('time, seconds')

plt.figure()
#plt.subplot(211)
plt.ylabel(r'degrees, degrees/second')
plt.plot(data['t'], data['qd0']*180./np.pi, label='$\dot{\psi}$')
plt.plot(data['t'], data['q0']*180./np.pi, label='$\psi$')
plt.legend(loc=0)
plt.xlabel('time, seconds')


plt.figure()
plt.plot(data['t'], data['u0'], label='$\omega_x$')
plt.plot(data['t'], data['u1'], label='$\omega_y$')
plt.plot(data['t'], data['u2'], label='$\omega_z$')
plt.title('Body fixed angular velocity')
plt.legend(loc=0)

plt.figure()
plt.title('Mechanical energy')
plt.subplot(311)
plt.plot(data['t'], data['ke'], label='ke')
plt.legend(loc=0)
plt.subplot(312)
plt.plot(data['t'], data['pe'], label='pe')
plt.legend(loc=0)
plt.subplot(313)
plt.plot(data['t'], data['te'] - data['te'][0], label=r'$\Delta E$')
plt.legend(loc=0)

plt.figure()
plt.subplot(211)
plt.plot(data['t'], data['Rx'], label='$\mathbf{F}_x$')
plt.legend(loc=0)
plt.plot(data['t'], data['Ry'], label='$\mathbf{F}_y$')
plt.legend(loc=0)
plt.subplot(212)
plt.plot(data['t'], data['Rz'], label='$\mathbf{F}_z$')
plt.legend(loc=0)
plt.title('Contact point reaction force')
plt.legend(loc=0)

plt.figure()
plt.plot(data['t'], data['ud0'], label='$\dot{\omega}_x$')
plt.plot(data['t'], data['ud1'], label='$\dot{\omega}_y$')
plt.plot(data['t'], data['ud2'], label='$\dot{\omega}_z$')
plt.title('Body fixed angular acceleration')
plt.legend(loc=0)

plt.figure()
plt.plot(data['q3'],data['q4'])
plt.title('Contact point location')

plt.show()

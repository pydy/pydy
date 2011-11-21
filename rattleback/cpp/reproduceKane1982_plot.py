import numpy as np
import matplotlib.pyplot as plt
import os
from math import pi

# Define a data type for the simdata structure
simdata = np.dtype([('t', np.float64),
                    ('alpha', np.float64),
                    ('beta', np.float64),
                    ('gamma', np.float64),
                    ('w1', np.float64),
                    ('w2', np.float64),
                    ('w3', np.float64),
                    ('delta', np.float64),
                    ('alpha1', np.float64),
                    ('alpha2', np.float64),
                    ('alpha3', np.float64),
                    ('ke', np.float64)])


os.system('rm -rf ./datafile.dat')
os.system('make')
os.system('./reproduceKane1982_sim')                # run the simulation
data = np.fromfile('datafile.dat', dtype=simdata) # read the data

plt.figure()
plt.subplot(211)
plt.plot(data['t'], data['delta']*180./pi, label=r'$\delta$')
plt.subplot(212)
plt.plot(data['t'], data['gamma']*180./pi, label=r'$\gamma$')
plt.title('Orientation')
plt.legend(loc=0)

plt.figure()
plt.plot(data['t'], data['w1'], label=r'$\omega_1$')
plt.plot(data['t'], data['w2'], label=r'$\omega_2$')
plt.plot(data['t'], data['w3'], label=r'$\omega_3$')
plt.title('Angular velocity')
plt.legend(loc=0)

plt.figure()
plt.plot(data['t'], data['alpha1'], label=r'$\alpha_1$')
plt.plot(data['t'], data['alpha2'], label=r'$\alpha_2$')
plt.plot(data['t'], data['alpha3'], label=r'$\alpha_3$')
plt.title('Angular acceleration')
plt.legend(loc=0)

plt.figure()
plt.plot(data['t'], data['ke'], label=r'ke')
plt.title('Kinetic energy')
plt.legend(loc=0)
plt.show()

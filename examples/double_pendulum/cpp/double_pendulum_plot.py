import numpy as np
import matplotlib.pyplot as plt
import os

# Define a data type for the simdata structure
simdata = np.dtype([('t', np.float64),
                    ('q1', np.float64),
                    ('q2', np.float64),
                    ('u1', np.float64),
                    ('u2', np.float64),
                    ('pe', np.float64),
                    ('ke', np.float64)])

os.system('./double_pendulum_sim')                # run the simulation
data = np.fromfile('datafile.dat', dtype=simdata) # read the data

plt.plot(data['t'], data['q1'], label='$q_1$')
plt.plot(data['t'], data['q2'], label='$q_2$')
plt.plot(data['t'], data['u1'], label='$u_1$')
plt.plot(data['t'], data['u2'], label='$u_2$')
plt.plot(data['t'], data['ke'], label='$ke$')
plt.plot(data['t'], data['pe'], label='$pe$')
plt.plot(data['t'], data['pe'] + data['ke'], label='$ke+pe$')
plt.legend(loc=0)

plt.show()

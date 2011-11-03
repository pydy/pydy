import numpy as np
import matplotlib.pyplot as plt
import os

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
                    ('Rx', np.float64),
                    ('Ry', np.float64),
                    ('Rz', np.float64),
                    ('ke', np.float64),
                    ('pe', np.float64),
                    ('te', np.float64),
                    ('delta', np.float64)])


os.system('rm -rf ./datafile.dat')
os.system('make')
os.system('./rattleback_sim')                # run the simulation
data = np.fromfile('datafile.dat', dtype=simdata) # read the data

plt.figure()
plt.subplot(211)
plt.plot(data['t'], data['q0'], label=r'$\alpha$')
plt.legend(loc=0)
plt.subplot(212)
plt.plot(data['t'], data['delta'], label=r'$\delta$')
plt.legend(loc=0)

plt.figure()
plt.plot(data['q3'],data['q4'])
plt.title('Contact point location')

plt.figure()
plt.title('Mechanical energy')
plt.subplot(311)
plt.plot(data['t'], data['ke'], label='ke')
plt.legend(loc=0)
plt.subplot(312)
plt.plot(data['t'], data['pe'], label='pe')
plt.legend(loc=0)
plt.subplot(313)
plt.plot(data['t'], data['te'], label='ke + pe')
plt.legend(loc=0)

#plt.plot(data['t'], data['q1'], label='Roll')
#plt.plot(data['t'], data['q2'], label='Pitch')
#plt.plot(data['t'], data['u0'], label='$u_0$')
#plt.plot(data['t'], data['u1'], label='$u_1$')
#plt.plot(data['t'], data['u2'], label='$u_2$')
plt.figure()
plt.subplot(211)
plt.plot(data['t'], data['Rx'], label='$\mathbf{F} \cdot \mathbf{y}_x$')
plt.legend(loc=0)
plt.plot(data['t'], data['Ry'], label='$\mathbf{F} \cdot \mathbf{y}_y$')
plt.legend(loc=0)
plt.subplot(212)
plt.plot(data['t'], data['Rz'], label='$\mathbf{F} \cdot \mathbf{y}_z$')
plt.legend(loc=0)
plt.title('Contact point reaction force')
plt.legend(loc=0)

plt.show()

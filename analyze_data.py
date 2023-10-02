import numpy as np
from matplotlib import pyplot as plt
import sys

import h5py

U_file = sys.argv[1]
h5f1 = h5py.File(U_file,'r')
mus = h5f1['mus'][:]

n_mus = mus.size

Gl1 = h5f1['Gl1_r'][:]

n_orbs = Gl1.shape[0]//2

n_particles = np.zeros(n_mus)
Delta_A = np.zeros(n_mus)
Delta_B = np.zeros(n_mus)
particles_A = np.zeros(n_mus)
particles_B = np.zeros(n_mus)

for n in range(n_mus):
    Gl = h5f1['Gl'+str(n+1)+'_r'][:]+1.0j*h5f1['Gl'+str(n+1)+'_i'][:]
    Deltas = np.zeros(n_orbs)
    particles = np.zeros(n_orbs)
    for m in range(n_orbs):
        n_particles[n] += 2.0*np.real(-1.0j*Gl[2*m,2*m])
        Deltas[m] = np.abs(-1.0j*Gl[2*m,2*m+1])
        particles[m] = np.real(-1.0j*Gl[2*m,2*m])
    particles_A[n] = particles[0]/2
    particles_B[n] = particles[1]/2
    Delta_A[n] = np.abs(Deltas[0])
    Delta_B[n] = np.abs(Deltas[1])

fig,ax = plt.subplots()
ax.plot(mus,n_particles)
fig,ax = plt.subplots()
ax.plot(mus,Delta_A)
ax.plot(mus,Delta_B)
fig,ax = plt.subplots()
ax.plot(mus,particles_A)
ax.plot(mus,particles_B)


plt.show()
